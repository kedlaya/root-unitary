/*
  Low-level code to exhaust over trees of Weil polynomials.
  This code does not implement parallelism; see the Cython wrapper.

#*****************************************************************************
#       Copyright (C) 2019-2026 Kiran S. Kedlaya <kskedl@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

*/

#include "power_sums.h"

/* Check for OpenMP at runtime. */
int num_threads() {
  #if defined(_OPENMP)
  return omp_get_max_threads();
  #endif
  return(1);
}

/*****
  Arithmetic functions

  As with FLINT library functions, aliasing is allowed unless specified.
*****/

inline int is_mpz(fmpz f) {
  return(COEFF_IS_MPZ(f));
}

/* Set res to floor((a+b)/2). */
inline void fmpz_fmid(fmpz_t res, const fmpz_t a, const fmpz_t b) {
  fmpz_add(res, a, b);
  fmpz_fdiv_q_ui(res, res, 2);
}

/* Set res to ceil((a+b)/2). */
inline void fmpz_cmid(fmpz_t res, const fmpz_t a, const fmpz_t b) {
  fmpz_add(res, a, b);
  fmpz_cdiv_q_ui(res, res, 2);
}

/* Set res to ceil(sqrt(a)). For the floor, use FLINT's built-in fmpz_sqrt instead. */
inline void fmpz_sqrt_c(fmpz_t res, const fmpz_t a) {
  int s = fmpz_root(res, a, 2);
  if (!s) fmpz_add_ui(res, res, 1);
}

/* Set res to floor((a/b + c sqrt(q))/d). No aliasing allowed. b and d must be positive.
   If b is NULL we interpret it as 1. */
inline void fmpq_floor_quad(fmpz_t res, const fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_t d, const fmpz_t q) {
  fmpz_mul(res, c, c);
  fmpz_mul(res, res, q);
  if (fmpz_sgn(c) >= 0) fmpz_sqrt(res, res);
  else {
    fmpz_sqrt_c(res, res);
    fmpz_neg(res, res);
  }
  if (b != NULL) fmpz_mul(res, res, b);
  fmpz_add(res, res, a);
  if (b != NULL) fmpz_fdiv_q(res, res, b);
  fmpz_fdiv_q(res, res, d);
}

/* Set res to ceil((a/b + c sqrt(q))/d). No aliasing allowed. b and d must be positive.
   If b is NULL we interpret it as 1. */
inline void fmpq_ceil_quad(fmpz_t res, const fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_t d, const fmpz_t q) {
  fmpz_mul(res, c, c);
  fmpz_mul(res, res, q);
  if (fmpz_sgn(c) >= 0) fmpz_sqrt_c(res, res);
  else {
    fmpz_sqrt(res, res);
    fmpz_neg(res, res);
  }
  if (b != NULL) fmpz_mul(res, res, b);
  fmpz_add(res, res, a);
  if (b != NULL) fmpz_cdiv_q(res, res, b);
  fmpz_cdiv_q(res, res, d);
}

/*
  Compute the Hankel determinant associated to the sequence seq of odd length n
  by a direct application of FLINT's determinant function.
*/

void hankel_determinant_direct(fmpz_t res, const fmpz *seq, int n) {
  fmpz_mat_t mat;
  int i, j, s;

  s = n/2 + 1;
  fmpz_mat_init(mat, s, s);
  for (i=0; i<s; i++)
    for (j=0; j<s; j++)
      fmpz_set(fmpz_mat_entry(mat, i, j), seq+i+j);
  fmpz_mat_det(res, mat);
  fmpz_mat_clear(mat);
  return;
}

/*
  Attempt to compute the Hankel determinant associated to the sequence seq of odd length n
  using Dodgson condensation. Returns 1 for success, 0 for failure.

  This function assumes that {w, 2*n-4} is scratch space.
*/

int hankel_determinant_condensation(fmpz_t res, const fmpz *seq, int n, fmpz *w) {
  if (n == 1) {
    fmpz_set(res, seq);
    return(1);
  }

  int i, n1 = n-2;
  fmpz *f0 = w; // Length n-2
  fmpz *f1 = w+n-2; // Length n-2
  fmpz *t = seq;

  for (i=0; i<n-2; i++) fmpz_fmms(f1+i, seq+i, seq+i+2, seq+i+1, seq+i+1);
  while (n1 > 1) {
    for (i=0; i<n1-2; i++) {
      if (fmpz_is_zero(t+i+2)) return(0); // Failure because of zero division
      fmpz_fmms(f0+i, f1+i, f1+i+2, f1+i+1, f1+i+1);
      fmpz_divexact(f0+i, f0+i, t+i+2);
    }
    n1 -= 2;
    t = f1; f1 = f0; f0 = t;
  }
  fmpz_set(res, f1);
  return(1);
}

/*
    Use a subresultant (Sturm-Habicht) sequence to test whether a given
    polynomial has all real roots. Note that this test has an early abort
    mechanism: having all real roots means that the sign sequence has
    the maximal number of sign changes, so the test aborts as soon
    as a sign change is missed.

    This function assumes that:
        - {poly, n} is a normalized vector with n >= 2 and positive leading coefficient;
        - {f0, n-1} and {f1, n-2} are scratch space.

    We add a (if b is NULL) or a*b (otherwise) to the constant term before testing.

    Based on code by Sebastian Pancratz from the FLINT repository (plus the Ducos variation).
*/

int _fmpz_poly_all_real_roots(fmpz *poly, long n, fmpz *f0, fmpz *f1, int force_squarefree,
			      const fmpz_t a, const fmpz_t b) {
  if (n <= 2) return(1);  // Constant or linear polynomial
  int sgn, i;

  /* Set f0 := deriv(poly). */
  _fmpz_poly_derivative(f0, poly, n);
  n -= 2;

  fmpz *t = f0+n;
  fmpz *t1 = poly+n+1;
  fmpz *content = NULL; // No content to remove in the first iteration

  /* Update the constant term of poly, result in f1. */
  if (b == NULL) fmpz_add(f1, a, poly); // Treat b as 1
  else {
    fmpz_mul(f1, a, b);
    fmpz_add(f1, f1, poly);
  }
  fmpz_mul(f1, f1, t);

  /* At this point deg(poly) = n+1, deg(f0) = n.
     We compute the leading coefficient of the pseudoremainder of poly modulo f0. */

  fmpz_fmms(f1+n, poly+n, t, f0+n-1, t1);
  if (n>1) fmpz_fmms(f1+n-1, poly+n-1, t, f0+n-2, t1);
  fmpz_fmms(f1+n-1, f1+n-1, t, f0+n-1, f1+n);

  /* If we miss any one sign change, we cannot have enough. */
  sgn = fmpz_sgn(f1+n-1);
  if (sgn == 1 || (force_squarefree && sgn == 0)) return(0);
  
  /* If f0 is a scalar, it is nonzero and we win. */
  if (n == 1) return(1);

  /* Set f1 to the pseudoremainder of poly modulo f0 in two steps:
       f1 --> f0[n+1]*poly - poly[n+2]*x*f0
       f1 --> f0[n+1]*f1 - f1[n+1]*f0.
     There is no Ducos variation at this step.
  */

  for (i=0; i<n-2; i++) fmpz_fmms(f1+i+1, poly+i+1, t, f0+i, t1);
  t1 = f1+n;
  for (i=0; i<n-1; i++) fmpz_fmms(f1+i, f1+i, t, f0+i, t1);

  while (1) {
    n--;

    /* We have deg(f0) == n+1, deg(f1) == n.
       If not forcing squarefree but sgn == 0, we win iff f0 = 0. */
    if (!force_squarefree && !sgn) return (_fmpz_vec_is_zero(f1, n));

    /* Divide f1 by content. */
    if (content != NULL) _fmpz_vec_scalar_divexact_fmpz(f1, f1, n+1, content);

    /* Compute the next leading coefficient times f0[n+1] and test the sign change. */
    t = f0+n-1;
    fmpz_fmms(t, t, f1+n, f1+n-1, f0+n);
    fmpz_divexact(t, t, f0+n+1);
    if (n>1) fmpz_sub(t, t, f1+n-2);
    fmpz_fmma(t, t, f1+n, f1+n-1, f1+n-1);
    sgn = fmpz_sgn(t);

    /* If we miss any one sign change, we cannot have enough. */
    if (sgn == 1 || (force_squarefree && sgn == 0)) return(0);

    /* If f0 is a scalar, it is nonzero and we win. */
    if (n == 1) return(1);

    /* We compute the pseudoremainder of f0 modulo f1 following a recipe of Ducos,
       leaving f0[n] and f0[n+1] intact.
       We also leave f0[n-1] intact as it has already been updated. */

    /* Take the remainder of f1[n]*(f0 (mod x^n)) modulo f1. */
    t = f1+n;
    t1 = f0+n;
    for (i=0; i<n-1; i++) fmpz_fmms(f0+i, f0+i, t, f1+i, t1);

    /* Divide (f0 mod x^{n-1}) by content. */
    content = f0+n+1;
    _fmpz_vec_scalar_divexact_fmpz(f0, f0, n-1, content);

    /* Subtract x*(f1 mod x^{n-2}) from f0. */
    _fmpz_vec_sub(f0+1, f0+1, f1, n-2);

    /* Multiply f0 by f1[n] and add f1[n-1]*(f1 mod x^{n-1}). */
    t1 = f1+n-1;
    for (i=0; i<n-1; i++) fmpz_fmma(f0+i, f0+i, t, f1+i, t1);

    /* Swap f0 with f1 at the pointer level. */
    t = f0; f0 = f1; f1 = t;
  }
}

/*****
  Memory allocation and deallocation
*****/

/* Static memory allocation and initialization. */
ps_static_data_t *ps_static_init(int d, fmpz_t q, fmpz_t lead, fmpz *modlist, 
                                 long node_limit, int force_squarefree) {
  int i, j, q_is_1, q_is_square;
  fmpz *k0, *pol;

  ps_static_data_t *st_data = (ps_static_data_t *)malloc(sizeof(ps_static_data_t));

  st_data->d = d;
  fmpz_init_set(st_data->q, q);
  q_is_1 = fmpz_is_one(q);
  q_is_square = fmpz_is_square(q);
  st_data->node_limit = node_limit;
  st_data->force_squarefree = force_squarefree;

  st_data->modlist = _fmpz_vec_init(d+1);
  k0 = st_data->modlist; // Used as a temporary variable for now

  st_data->sum_mats = _fmpz_vec_init((d+1)*(d+1));
  for (i=0; i<=d; i++) {
    /* Coefficients of 2*(i-th Chebyshev polynomial)(x/2).
       The coefficient of x^j is multiplied by c^{i-j} q^{(i-j)/2} where c is the leading coefficient. */
    pol = st_data->sum_mats + (d+1)*i;
    _fmpz_poly_chebyshev_t(pol, i);
    _fmpz_poly_scale_2exp(pol, i+1, -1);
    if (!q_is_1)
      for (j=i%2; j<=i; j+=2) {
        fmpz_pow_ui(k0, q, (i-j)/2);
        fmpz_mul(pol+j, k0, pol+j);
      }
    if (!fmpz_is_one(lead))
      for (j=i%2; j<=i; j+=2) {
        fmpz_pow_ui(k0, lead, i-j);
        fmpz_mul(pol+j, k0, pol+j);
      }
  }

  for (i=0; i<=d; i++) fmpz_set(st_data->modlist + i, modlist+d-i);

  st_data->binom_mat = _fmpz_vec_init((d+1)*(d+1));
  for (i=0; i<=d; i++)
    for (j=0; j<=d; j++)
      fmpz_bin_uiui(st_data->binom_mat+(d+1)*i+j, j, i);

  st_data->eval_pm2_mats = _fmpz_vec_init(2*(d+1)*(d+1));
  for (i=0; i<=d; i++)
    for (j=0; j<=i; j++) {
      k0 = st_data->eval_pm2_mats + (d+1)*(2*i+j%2) + j;
      if (q_is_1) fmpz_one(k0);
      else if (q_is_square) {
        fmpz_sqrt(k0, q);
        fmpz_pow_ui(k0, k0, j);
      } else fmpz_pow_ui(k0, q, j/2);
      fmpz_mul_2exp(k0, k0, j);
      fmpz_mul_si(k0, k0, -i);
      fmpz_mul(k0, k0, st_data->binom_mat + (d+2)*(d-i) + j);
    }

  return(st_data);
}

/* Dynamic memory allocation and initialization.
   Call with coefflist == NULL to prepare an inactive process. */
ps_dynamic_data_t *ps_dynamic_init(int d, fmpz_t q, fmpz *coefflist) {
  ps_dynamic_data_t *dy_data = (ps_dynamic_data_t *)malloc(sizeof(ps_dynamic_data_t));

  dy_data->d = d;
  dy_data->n = d;
  dy_data->node_count = 0;
  dy_data->ascend = 0;
  dy_data->pol = _fmpz_vec_init(d+1);
  dy_data->sympol = _fmpz_vec_init(2*d+1);
  if (coefflist != NULL) {
    dy_data->flag = 1; // Activate this process
    _fmpz_vec_set(dy_data->pol, coefflist, d+1);
  } else dy_data->flag = 0; // Flag this process as inactive
  dy_data->upper = _fmpz_vec_init(d+1);
  dy_data->power_sums_num = _fmpz_vec_init(d+1);
  fmpz_set_si(dy_data->power_sums_num, d);
  _fmpz_vec_zero(dy_data->power_sums_num+1, d);

  dy_data->hankel_dets = _fmpz_vec_init(2*d+2);
  fmpz_set_si(dy_data->hankel_dets, d);
  fmpz_set_si(dy_data->hankel_dets+1, 1);

  dy_data->wlen = 3*d+10;
  dy_data->w = _fmpz_vec_init(dy_data->wlen);

  return(dy_data);
}

/* Static memory deallocation. */
void ps_static_clear(ps_static_data_t *st_data) {
  if (st_data == NULL) return;

  int d = st_data->d;

  fmpz_clear(st_data->q);
  _fmpz_vec_clear(st_data->modlist, d+1);
  _fmpz_vec_clear(st_data->binom_mat, (d+1)*(d+1));
  _fmpz_vec_clear(st_data->sum_mats, (d+1)*(d+1));
  _fmpz_vec_clear(st_data->eval_pm2_mats, 2*(d+1)*(d+1));
  free(st_data);
}

/* Dynamic memory deallocation. */
void ps_dynamic_clear(ps_dynamic_data_t *dy_data) {
  if (dy_data == NULL) return;

  int d = dy_data->d;

  _fmpz_vec_clear(dy_data->pol, d+1);
  _fmpz_vec_clear(dy_data->sympol, 2*d+1);
  _fmpz_vec_clear(dy_data->upper, d+1);
  _fmpz_vec_clear(dy_data->power_sums_num, d+1);
  _fmpz_vec_clear(dy_data->hankel_dets, 2*d+2);
  _fmpz_vec_clear(dy_data->w, dy_data->wlen);
  free(dy_data);
}

/*****
  Low-level flow control
*****/

/* Increment the current moving counter and update stored data to match. 
   If step is NULL it is interpreted as 1. */

inline void step_forward(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int n, fmpz_t step) {
  int d = st_data->d;
  int k = d - n;
  fmpz *pol = dy_data->pol;
  fmpz *poln = pol+n;
  fmpz *modulus = st_data->modlist + n;
  int modulus_is_1 = fmpz_is_one(modulus);
  fmpz *pow_num = dy_data->power_sums_num+k;
  fmpz *det = dy_data->hankel_dets + 2*k;
  fmpz *t0z = dy_data->w;

  if (fmpz_is_one(pol+d)) fmpz_set_ui(t0z, k);
  else {
    fmpz_pow_ui(t0z, pol+d, k-1);
    fmpz_mul_ui(t0z, t0z, k);
  }
  if (!modulus_is_1) fmpz_mul(t0z, t0z, modulus);
  if (step == NULL)
    if (modulus_is_1) fmpz_add_ui(poln, poln, 1);
    else fmpz_add(poln, poln, modulus);
  else {
    fmpz_mul(t0z, t0z, step);
    if (modulus_is_1) fmpz_add(poln, poln, step);
    else fmpz_addmul(poln, step, modulus);
  }
  fmpz_sub(pow_num, pow_num, t0z);
  if (k > 1) {
    fmpz_submul(det, det-4, t0z);
    fmpz_addmul(det+1, det-3, t0z);
  } else {
    fmpz_sub(det, det, t0z);
    fmpz_add(det+1, det+1, t0z);
  }
}

/* Given a polynomial pol of length k, shrink the interval [lower, upper] to the range of constant terms
   which when added to pol give a polynomial with all real roots. Returns 0 if this range is empty (with
   lower, upper now undefined) and 1 otherwise (with lower, upper now the endpoints of the new range).

   This function assumes that {w, 2*k+2} is scratch space.
*/

int apply_rolle_condition(fmpz_t lower, fmpz_t upper, const fmpz *pol, int k, int force_squarefree, const fmpz_t modulus, fmpz *w) {
  int r, s;

  /* Allocate from working space. */
  fmpz *t0z = w;
  fmpz *t1z = w+1;
  fmpz *t2z = w+2;
  fmpz *f0 = w+3; // Length k
  fmpz *f1 = w+k+3; // Length k-1

  #define TEST_ROOTS(x) _fmpz_poly_all_real_roots(pol, k, f0, f1, force_squarefree, x, modulus)

  /* Handle the case upper == lower directly. */
  fmpz_sub(t0z, upper, lower);
  if (fmpz_is_zero(t0z)) return(TEST_ROOTS(lower));

  /* Look for a single value where the Rolle criterion holds. */
  fmpz_add_ui(t0z, t0z, 1);
  r = fmpz_flog_ui(t0z, 2); // r = floor(log_2 (upper-lower+1)); forced to be positive
  do {
    fmpz_one_2exp(t2z, r);
    if (!r) fmpz_set(t0z, lower);
    else {
      fmpz_add(t0z, lower, t2z);
      fmpz_sub_ui(t0z, t0z, 1);
    }
    do {
      if (s = TEST_ROOTS(t0z)) break;
      else fmpz_addmul_ui(t0z, t2z, 2);
    } while (fmpz_cmp(t0z, upper) <= 0);
    if (s) break;
    else r--;
  } while (r >= 0);
  if (!s) return(0); // Found nothing

  if (r == 0) { // In this case, enforce lower == upper and exit
    fmpz_set(lower, t0z);
    fmpz_set(upper, t0z);
    return(1);
  }

  /* Shorten the interval based on tested values. */
  fmpz_sub(t1z, t0z, t2z);
  fmpz_add_ui(lower, t1z, 1); // Does not decrease lower
  fmpz_add(t1z, t0z, t2z);
  if (fmpz_cmp(t1z, upper) <= 0) fmpz_sub_ui(upper, t1z, 1);

  /* Use binary searches to compute the interval on which the Rolle criterion is satisfied. */
  fmpz_set(t1z, t0z);
  while (fmpz_cmp(lower, t0z)) {
    fmpz_fmid(t2z, lower, t0z);
    if (TEST_ROOTS(t2z)) fmpz_set(t0z, t2z);
    else fmpz_add_ui(lower, t2z, 1);
  }
  while (fmpz_cmp(t1z, upper)) {
    fmpz_cmid(t0z, t1z, upper);
    if (TEST_ROOTS(t0z)) fmpz_set(t1z, t0z);
    else fmpz_sub_ui(upper, t0z, 1);
  }
  return(1);
}

/* The following is the key subroutine: given some initial coefficients, compute
   a lower and upper bound for the next coefficient. Return 1 iff the resulting
   interval is nonempty.

   The value of dy_data->pol+n is assumed to be correct modulo st_data->modlist+n.
*/
int set_range_from_power_sums(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int n) {
  if (n < 0) return(1);

  /* Static data */

  int d = st_data->d;
  int k = d - n;
  int force_squarefree = st_data->force_squarefree;
  fmpz *modulus = st_data->modlist + n;
  int modulus_is_0 = fmpz_is_zero(modulus);
  int modulus_is_1 = fmpz_is_one(modulus);
  fmpz *q = st_data->q;
  int q_is_1 = fmpz_is_one(q);
  int q_is_square = q_is_1 ? 1 : fmpz_is_square(q);

  /* Dynamic data */

  fmpz *pol = dy_data->pol + n;
  fmpz *lead = pol + k;
  int lead_is_1 = fmpz_is_one(lead);
  fmpz *pow_num = dy_data->power_sums_num;

  /* Integers allocated from working space, maintained throughout */

  fmpz *f = dy_data->w;
  fmpz *lead_pow = f+1;
  fmpz_pow_ui(lead_pow, lead, k-1);

  fmpz *upper = f+2;
  fmpz *lower = f+3;

  /* Local working space, not maintained throughout */

  int i, s;
  fmpz *w = f+4; // Length 3*k+3, passed to subroutines

  fmpz *t0z = w;
  fmpz *t1z = w+1;
  fmpz *t2z = w+2; // Affected by change_by_sign

  /* Unallocated pointers */

  fmpz *tz, *tza;

  /* Adjust lower and upper bounds within set_range_from_power_sums (affects lower, upper, t2z).
     The pair (val1, val2) stands for g = val1 + val2*sqrt(q).
     The value in val1 is specified as a numerator-denominator
     pair which need not be canonicalized (a denominator of NULL is interpreted as 1).
     A value of NULL for val2_num is interpreted as 0.
     No aliasing allowed unless val2_num = NULL, in which case allowed between t2z and val1_num.

     Given that g is a monic linear function of the k-th power sum, then:

     -- passing r = 0 imposes the condition g >= 0 (or g > 0 if force_squarefree != 0);
     -- passing r = 1 imposes the condition g <= 0 (or g < 0 if force_squarefree != 0);
     -- passing update = 0 means we are setting the initial bounds;
     -- passing update = 1 means we are updating previously set bounds.
  */

  inline void change_by_sign(int update, int r, const fmpz_t val1_num, const fmpz_t val1_den, const fmpz_t val2_num) {
    fmpz *k;

    if (val2_num == NULL) {
      if (r ^ force_squarefree) {
        fmpz_cdiv_q(t2z, val1_num, f);
        if (val1_den != NULL) fmpz_cdiv_q(t2z, t2z, val1_den);
      } else {
        fmpz_fdiv_q(t2z, val1_num, f);
        if (val1_den != NULL) fmpz_fdiv_q(t2z, t2z, val1_den);
      }
    } else if (r ^ force_squarefree) fmpq_ceil_quad(t2z, val1_num, val1_den, val2_num, f, q);
    else fmpq_floor_quad(t2z, val1_num, val1_den, val2_num, f, q);
    if (!r) { // change_upper
      if (force_squarefree) {
        if (!update || fmpz_cmp(t2z, upper) <= 0) fmpz_sub_ui(upper, t2z, 1);
      } else if (!update || fmpz_cmp(t2z, upper) < 0) fmpz_set(upper, t2z);
    } else { // change_lower
      if (force_squarefree) {
        if (!update || fmpz_cmp(t2z, lower) >= 0) fmpz_add_ui(lower, t2z, 1);
      } else if (!update || fmpz_cmp(t2z, lower) > 0) fmpz_set(lower, t2z);
    }
  }

  /* If modulus==0, reduce the interval to [0]. */

  fmpz_set_ui(f, k);
  if (modulus_is_0) {
    fmpz_zero(lower);
    fmpz_zero(upper);
  } else if (!modulus_is_1) fmpz_mul(f, f, modulus);

  /* Update pow_num[k] using the Girard-Newton formula.
     This is the k-th power sum times c^k where c is the leading coefficient. */

  fmpz_set_ui(pow_num, k); // Temporary change to apply Girard-Newton
  fmpz_zero(pow_num+k);
  for (i=0; i<k; i++) {
    if (!lead_is_1 && i > 0) fmpz_mul(t0z, t1z, pol+k-1-i);
    else fmpz_set(t0z, pol+k-i-1);
    fmpz_submul(pow_num+k, t0z, pow_num+k-1-i);
    if (!lead_is_1)
      if (i == 0) fmpz_set(t1z, lead);
      else if (i < k-1) fmpz_mul(t1z, t1z, lead);
  }
  fmpz_set_ui(pow_num, d); // Change back to the correct value, needed for Chebyshev criterion

  /* Chebyshev criterion: the k-th symmetrized power sum must lie in [-2*d*q^(k/2), 2*d*q^(k/2)]. */

  if (lead_is_1) fmpz_set_ui(t1z, 2*d);
  else fmpz_mul_ui(t1z, lead, 2*d);
  if (!q_is_1 && k > 1) {
    fmpz_pow_ui(t0z, q, k/2);
    fmpz_mul(t1z, t1z, t0z);
  }
  if (!q_is_1 && q_is_square && k%2 == 1) {
    fmpz_sqrt(t0z, q);
    fmpz_mul(t1z, t1z, t0z);
  }
  _fmpz_vec_dot(t0z, st_data->sum_mats+(d+1)*k, pow_num, k+1);
  if (q_is_square || k%2 == 0) { // q^{k/2} is rational
    if (!lead_is_1) fmpz_mul(t1z, lead_pow, t1z);
    fmpz_add(t2z, t0z, t1z);
    change_by_sign(modulus_is_0, 0, t2z, lead_pow, NULL);
    fmpz_sub(t2z, t0z, t1z);
    change_by_sign(modulus_is_0, 1, t2z, lead_pow, NULL);
  } else { // q^{k/2} is irrational
    change_by_sign(modulus_is_0, 0, t0z, lead_pow, t1z);
    fmpz_neg(t1z, t1z);
    change_by_sign(modulus_is_0, 1, t0z, lead_pow, t1z);
  }

  /* Descartes criterion: the evaluations of the n-th derivative of pol at -2*sqrt(q), 2*sqrt(q)
     have the correct signs. */

  _fmpz_vec_dot(t0z, st_data->eval_pm2_mats+(d+1)*(2*k), pol, k+1);
  _fmpz_vec_dot(t1z, st_data->eval_pm2_mats+(d+1)*(2*k+1), pol, k+1);
  if (q_is_square) {
    fmpz_add(t2z, t0z, t1z);
    change_by_sign(1, 1, t2z, NULL, NULL);
    fmpz_sub(t2z, t0z, t1z);
    change_by_sign(1, 1-(k%2), t2z, NULL, NULL);
  } else {
    change_by_sign(1, 1, t0z, NULL, t1z);
    fmpz_neg(t1z, t1z);
    change_by_sign(1, 1-(k%2), t0z, NULL, t1z);
  }
  if (fmpz_cmp(lower, upper) > 0) return(0);

  /* Hausdorff criterion: the relevant Hankel matrices have nonnegative determinant.
     In order to restrict to integer arithmetic, we skip this condition when
     k is odd and q is not a perfect square. */

  for (i=0; i<=1; i++) {
    if (!q_is_square && k%2 == 1) continue; // Hankel matrix is not defined over Q

    /* Build the sequence of entries of the appropriate Hankel matrix. */
    if (i == 1 && k%2 == 0) { // t0z := 4*lead^2*q
      s = k - 1;
      if (lead_is_1) fmpz_one(t0z);
      else fmpz_mul(t0z, lead, lead);
      if (!q_is_1) fmpz_mul(t0z, t0z, q);
      fmpz_mul_ui(t0z, t0z, 4);
    } else {
      s = k + 1 - k%2;
      if (k%2 == 1) // t0z := 2*lead*sqrt(q)
        if (q_is_1) fmpz_mul_ui(t0z, lead, 2);
        else {
          fmpz_sqrt(t0z, q);
          if (!lead_is_1) fmpz_mul(t0z, t0z, lead);
          fmpz_mul_ui(t0z, t0z, 2);
        }
    }
    if (i == 0 && k%2 == 0) tza = pow_num;
    else {
      tza = w;
      _fmpz_vec_scalar_mul_fmpz(tza, pow_num, s, t0z);
      if (k%2 == 0) _fmpz_vec_sub(tza, tza, pow_num+2, s);
      else if (i == 0) _fmpz_vec_add(tza, tza, pow_num+1, s);
      else _fmpz_vec_sub(tza, tza, pow_num+1, s);
    }

    /* Compute the determinant, then deduce a condition. */
    tz = dy_data->hankel_dets + 2*k + i;
    if (!hankel_determinant_condensation(tz, tza, s, w+k+1)) 
      hankel_determinant_direct(tz, tza, s);

    if (k > 1 && (force_squarefree || fmpz_sgn(tz-4) > 0)) {
      fmpz_set(t0z, tz);
      if (lead_is_1) fmpz_set(t1z, tz-4);
      else fmpz_mul(t1z, tz-4, lead_pow);
    } else { // If the determinant vanishes, argue that the corner entry is nonnegative
      fmpz_set(t0z, tza+s-1);
      if (lead_is_1) fmpz_one(t1z);
      else fmpz_set(t1z, lead_pow);
    }
    if (i == 1) fmpz_neg(t0z, t0z);
    change_by_sign(1, i, t0z, t1z, NULL);
  }
  if (fmpz_cmp(lower, upper) > 0) return(0);

  /* Rolle criterion: the divided n-th derivative of pol has all roots real. */

  tza = st_data->binom_mat + (d+2)*n;
  for (i=0; i<=k; i++) fmpz_mul(w+i, tza+i, pol+i);
  tz = modulus_is_1 ? NULL : modulus;
  if (!apply_rolle_condition(lower, upper, w, k+1, force_squarefree, tz, w+k+1)) return(0);

  /* Set the new upper bound. */

  if (!modulus_is_1) fmpz_mul(upper, upper, modulus);
  fmpz_add(dy_data->upper+n, pol, upper);

  /* Set the new polynomial value, then correct the k-th power sum and related quantities. */
  step_forward(st_data, dy_data, n, lower);

  return(1);
}

/*****
  High-level flow control
*****/

/* Split off a subtree (without allocating new memory).
   The donor process yields its current branch at the first coefficient that is not uniquely specified.
   The donee process may in turn be split immediately.

   It is safe to run this in parallel as long as the instances of dy_data are pairwise distinct
   and the instances of dy_data2 are pairwise distinct.
*/
void ps_dynamic_split(const ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, ps_dynamic_data_t *dy_data2) {
  if (dy_data == NULL || dy_data2 == NULL || dy_data->flag <= 0 || dy_data2->flag) return;

  int i;
  int d = dy_data->d;
  int n = dy_data->n;
  int ascend = dy_data->ascend;
  int k = n + ascend;

  for (i=d; i>k; i--)
    if (fmpz_cmp(dy_data->pol+i, dy_data->upper+i) < 0) {
      /* Copy the current state of the donor to the donee process. */
      int j = d - i + 1;
      _fmpz_vec_set(dy_data2->pol+i, dy_data->pol+i, j);
      _fmpz_vec_set(dy_data2->upper+i, dy_data->upper+i, j);
      _fmpz_vec_set(dy_data2->power_sums_num, dy_data->power_sums_num, j);
      _fmpz_vec_set(dy_data2->hankel_dets, dy_data->hankel_dets, 2*j);

      fmpz *lower = dy_data->pol+i;
      fmpz *upper = dy_data->upper+i;
      fmpz *modulus = st_data->modlist + i;
      int modulus_is_1 = fmpz_is_one(modulus);
      fmpz *t0z = dy_data->w; // Temporary variable

      /* Restrict the donee process to the right half of the interval. */
      fmpz_sub(t0z, upper, lower);
      if (!modulus_is_1) fmpz_divexact(t0z, t0z, modulus);
      fmpz_cdiv_q_ui(t0z, t0z, 2);
      step_forward(st_data, dy_data2, i, t0z);
      dy_data2->ascend = 0;
      dy_data2->n = i;
      dy_data2->flag = 1;

      /* Restrict the donor process to the left half of the interval. */
      fmpz_sub_ui(t0z, t0z, 1);
      if (!modulus_is_1) fmpz_mul(t0z, t0z, modulus);
      fmpz_add(upper, lower, t0z);
      break;
    }
}

/* Top-level flow control: allow one process to run for up to max_steps iterations,
   or until it finds a polynomial to be returned, whichever comes first.

   Return value also sent back in dy_data->flag:
   1: in process
   2: found a solution (returned in dy_data->sympol)
   0: tree exhausted
   -1: maximum number of nodes reached

   It is *not* threadsafe to run next_pol and dynamic_split simultaneously on the same dy_data.
*/

int next_pol(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int max_steps) {
  if (dy_data==NULL || !dy_data->flag) return(0); // No work assigned to this process
  dy_data->flag = 0; // Prevent work-stealing while this process is running

  int d = st_data->d;
  int n = dy_data->n;
  if (n>d) return(0); // This process is exhausted.

  int ascend = dy_data->ascend;
  long node_limit = st_data->node_limit;
  long node_count = dy_data->node_count;
  fmpz *pol = dy_data->pol;
  fmpz *upper = dy_data->upper;

  int i, j;
  int flag = 1;
  int count_steps = 0;

  while (count_steps <= max_steps) {
    count_steps += 1;
    if (ascend) { // Ascend the tree and step forward as needed.
      n += ascend;
      if (n > d) {
        flag = 0; // This process is exhausted.
        break;
      } else if (!(ascend = (fmpz_cmp(pol+n, upper+n) >= 0)))
	step_forward(st_data, dy_data, n, NULL);
    } else if (n < 0) { // Return a solution.
      fmpz *sympol = dy_data->sympol;
      fmpz *t;
      fmpz *q = st_data->q;
      int q_is_1 = fmpz_is_one(q);

      /* Convert pol into its reciprocal transform, stored in sympol. */
      for (i=0; i<=d; i++) {
        t = sympol + d - i;
        fmpz_set(t, pol+i);
	for (j=i; j>0; j--) {
	  if (j==i) fmpz_set(t+2*i, t);
	  else fmpz_add(t+2*j, t+2*j, t);
          if (!q_is_1) fmpz_mul(t, t, q);
	  fmpz_mul_ui(t, t, j);
	  fmpz_divexact_ui(t, t, i-j+1);
	}
      }
      ascend = 1;
      flag = 2;
      break;
    } else { // Compute children of the current node.
      n -= 1;
      if (ascend = !set_range_from_power_sums(st_data, dy_data, n)) { // Found a terminal node
	node_count += 1;
	if (node_limit != -1 && node_count >= node_limit) {
	  flag = -1;
	  break;
	}
      }
    }
  }
  /* Record the final working state. */
  dy_data->ascend = ascend;
  dy_data->n = n;
  dy_data->node_count = node_count;
  dy_data->flag = flag;
  return(flag);
}
