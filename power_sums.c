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

/* Check for OpenMP at runtime.
*/
int has_openmp() {
  #if defined(_OPENMP)
  return(1);
  #endif
  return(0);
}

int num_threads() {
  if (has_openmp()) return omp_get_max_threads(); 
  else return(1);
}

/*****
  Arithmetic functions

  As with FLINT library functions, aliasing is allowed unless specified.

  Unlike for FLINT library functions, input fmpq's need not be canonicalized, and
  output fmpq's are not guaranteed to be canonicalized.
*****/

inline int is_mpz(fmpz f) {
  return(COEFF_IS_MPZ(f));
}

/* Set res to a/b. No aliasing allowed. */
inline void fmpq_div_raw(fmpq_t res, const fmpq_t a, const fmpq_t b) {
  fmpz_mul(fmpq_numref(res), fmpq_numref(a), fmpq_denref(b));
  fmpz_mul(fmpq_denref(res), fmpq_denref(a), fmpq_numref(b));
}

/* Set res to floor(a). */
inline void fmpq_floor(fmpz_t res, const fmpq_t a) {
  fmpz_fdiv_q(res, fmpq_numref(a), fmpq_denref(a));
};

/* Set res to ceil(a). */
inline void fmpq_ceil(fmpz_t res, const fmpq_t a) {
  fmpz_cdiv_q(res, fmpq_numref(a), fmpq_denref(a));
};

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

/* Set res to floor(sqrt(a)). */
inline void fmpz_sqrt_f(fmpz_t res, const fmpz_t a) {
  fmpz_sqrt(res, a);
}

/* Set res to ceil(sqrt(a)). */
inline void fmpz_sqrt_c(fmpz_t res, const fmpz_t a) {
  int s = fmpz_root(res, a, 2);
  if (!s) fmpz_add_ui(res, res, 1);
}

/* Set res to floor(a + b sqrt(q)). No aliasing allowed. */
inline void fmpq_floor_quad(fmpz_t res, const fmpq_t a, const fmpq_t b, const fmpz_t q) {
  fmpz *anum = fmpq_numref(a);
  fmpz *aden = fmpq_denref(a);
  fmpz *bnum = fmpq_numref(b);
  fmpz *bden = fmpq_denref(b);
  int bden_s = fmpz_sgn(bden);
  int aden_s = fmpz_sgn(aden);
  int bnum_s = fmpz_sgn(bnum);

  fmpz_mul(res, aden, bnum);
  fmpz_mul(res, res, res);
  fmpz_mul(res, res, q);
  if (bnum_s*bden_s >= 0) fmpz_sqrt_f(res, res);
  else {
    fmpz_sqrt_c(res, res);
    fmpz_neg(res, res);
  }
  fmpz_mul_si(res, res, aden_s*bden_s);
  fmpz_addmul(res, anum, bden);
  if (bden_s > 0) fmpz_fdiv_q(res, res, aden);
  else fmpz_cdiv_q(res, res, aden);
  fmpz_fdiv_q(res, res, bden);
}

/* Set res to ceil(a + b sqrt(q)). No aliasing allowed. */
inline void fmpq_ceil_quad(fmpz_t res, const fmpq_t a, const fmpq_t b, const fmpz_t q) {
  fmpz *anum = fmpq_numref(a);
  fmpz *aden = fmpq_denref(a);
  fmpz *bnum = fmpq_numref(b);
  fmpz *bden = fmpq_denref(b);
  int bden_s = fmpz_sgn(bden);
  int aden_s = fmpz_sgn(aden);
  int bnum_s = fmpz_sgn(bnum);
    
  fmpz_mul(res, aden, bnum);
  fmpz_mul(res, res, res);
  fmpz_mul(res, res, q);
  if (bnum_s*bden_s >= 0) fmpz_sqrt_c(res, res);
  else {
    fmpz_sqrt_f(res, res);
    fmpz_neg(res, res);
  }
  fmpz_mul_si(res, res, aden_s*bden_s);
  fmpz_addmul(res, anum, bden);
  if (bden_s > 0) fmpz_cdiv_q(res, res, aden);
  else fmpz_fdiv_q(res, res, aden);
  fmpz_cdiv_q(res, res, bden);
}

/*
    Use a subresultant (Sturm-Habicht) sequence to test whether a given
    polynomial has all real roots. Note that this test has an early abort
    mechanism: having all real roots means that the sign sequence has
    the maximal number of sign changes, so the test aborts as soon
    as a sign change is missed.

    This function assumes that:
        - {poly, n} is a normalized vector with n >= 2 and nonzero leading coefficient
        - {f0, n} and {f1, n-1} are scratch space.
    If a and b are not NULL, we add a*b to the constant term before testing.

    Based on code by Sebastian Pancratz from the FLINT repository.
*/

int _fmpz_poly_all_real_roots(fmpz *poly, long n, fmpz *f0, fmpz *f1, int force_squarefree,
			      const fmpz_t a, const fmpz_t b) {
  if (n <= 2) return(1);  // Constant or linear polynomial
  fmpz *t;
  int i;

  /* Set f1 := deriv(poly). */
  _fmpz_poly_derivative(f1, poly, n);

  n--; // now n = deg(poly)
  int n0 = n; // Initial degree
  int sgn0_l = fmpz_sgn(poly+n); // Sign of initial leading coefficient

  while (1) { // At this point deg(f0) = n, deg(f1) = n-1.
    /* 
       We compute the pseudoremainder of f0 modulo f1 in two steps:
       f0 --> f1[n-1]*f0 - f0[n]*x*f1
       f0 --> f1[n-1]*f0 - f0[n-1]*f1
    */
    
    i = (n == n0);
    t = i ? poly : f0; // if n == n0, f0 is not yet initialized
    _fmpz_vec_scalar_mul_fmpz(f0+i, t+i, n-i, f1+n-1);
    _fmpz_vec_scalar_submul_fmpz(f0+1, f1, n-1, t+n);
    if (i) // Fill in constant term of f0 using a and b
      if (a != NULL && b != NULL) {
        fmpz_mul(f0, a, b);
        fmpz_fmma(f0, f0, f1+n-1, poly, f1+n-1);
      }
      else fmpz_mul(f0, poly, f1+n-1);

    n--; // At this point deg(f0) = deg(f1) = n.
    _fmpz_vec_scalar_mul_fmpz(f0, f0, n, f1+n);
    _fmpz_vec_scalar_submul_fmpz(f0, f1, n, f0+n);

    /* At this point deg(f0) = n-1, deg(f1) = n.
       If f0 = 0, we win unless we are insisting on squarefree. */
    if (!force_squarefree && _fmpz_vec_is_zero(f0, n)) return(1);

    /* If we miss any one sign change, we cannot have enough. */
    if ((n0-n) % 2 == 1) sgn0_l = -sgn0_l;
    if (fmpz_sgn(f0+n-1) != sgn0_l) return(0);

    /* If f0 is a scalar, it is nonzero and we win. */
    if (n == 1) return(1);

    /* Extract content from f0.
       This seems to do better in practice than an explicit subresultant computation.
       Note that f0+n is now allocated but unused, so available for a temporary value. */
    _fmpz_vec_content(f0+n, f0, n);
    _fmpz_vec_scalar_divexact_fmpz(f0, f0, n, f0+n);

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
  int i, j, q_is_1;
  ps_static_data_t *st_data;
  fmpz *k0, *pol;
  fmpq *k1;

  st_data = (ps_static_data_t *)malloc(sizeof(ps_static_data_t));

  st_data->d = d;
  fmpz_init_set(st_data->q, q);
  q_is_1 = fmpz_is_one(q);
  st_data->node_limit = node_limit;
  st_data->force_squarefree = force_squarefree;

  st_data->modlist = _fmpz_vec_init(d+1);
  st_data->f = _fmpq_vec_init(d+1);
  k0 = st_data->modlist; // Used as a temporary variable for now

  st_data->sum_mats = _fmpz_vec_init((d+1)*(d+1));
  for (i=0; i<=d; i++) {
    /* Coefficients of 2*(i-th Chebyshev polynomial)(x/2).
       If q != 1, the coeff of x^j is multiplied by q^{(i-j)/2}. */
    pol = st_data->sum_mats + (d+1)*i;
    _fmpz_poly_chebyshev_t(pol, i);
    _fmpz_poly_scale_2exp(pol, i+1, -1);
    if (!q_is_1)
      for (j=i%2; j<=i; j+=2) {
        fmpz_pow_ui(k0, q, (i-j)/2);
        fmpz_mul(pol+j, k0, pol+j);
      }
  }
 
  for (i=0; i<=d; i++) {
    k0 = st_data->modlist + i;
    k1 = st_data->f + i;
    fmpz_set(k0, modlist+d-i);
    fmpq_set_si(k1, d-i, 1);
    fmpq_div_fmpz(k1, k1, lead);
    /* In order to apply the Chebyshev and Descartes criteria
       when the modulus is 0, we must pretend that the modulus is 1. */
    if (!fmpz_is_zero(k0)) fmpq_mul_fmpz(k1, k1, k0);
  }

  st_data->binom_mat = _fmpz_vec_init((d+1)*(d+1));
  for (i=0; i<=d; i++)
    for (j=0; j<=d; j++)
      fmpz_bin_uiui(st_data->binom_mat+(d+1)*i+j, j, i);

  st_data->eval_pm2_mats = _fmpz_vec_init(2*(d+1)*(d+1));
  for (i=0; i<=d; i++)
    for (j=0; j<=i; j++) {
      k0 = st_data->eval_pm2_mats + (d+1)*(2*i+j%2) + j;
      if (q_is_1) fmpz_one(k0);
      else fmpz_pow_ui(k0, st_data->q, j/2);
      fmpz_mul_2exp(k0, k0, j);
      fmpz_mul_si(k0, k0, -i);
      fmpz_mul(k0, k0, st_data->binom_mat + (d+2)*(d-i) + j);
    }

  return(st_data);
}

/* Dynamic memory allocation and initialization.
   Call with coefflist == NULL to prepare an inactive process. */
ps_dynamic_data_t *ps_dynamic_init(int d, fmpz_t q, fmpz *coefflist) {
  ps_dynamic_data_t *dy_data;
  int i;

  dy_data = (ps_dynamic_data_t *)malloc(sizeof(ps_dynamic_data_t));
  dy_data->d = d;
  dy_data->n = d;
  dy_data->node_count = 0;
  dy_data->ascend = 0;
  dy_data->pol = _fmpz_vec_init(d+1);
  dy_data->sympol = _fmpz_vec_init(2*d+3);
  if (coefflist != NULL) {
    dy_data->flag = 1; // Activate this process
    for (i=0; i<=d; i++)
      fmpz_set(dy_data->pol+i, coefflist+i);
  } else dy_data->flag = 0; // Flag this process as inactive
  dy_data->upper = _fmpz_vec_init(d+1);

  fmpq_mat_init(dy_data->power_sums, d+1, 1);
  fmpq_set_si(fmpq_mat_entry(dy_data->power_sums, 0, 0), d, 1);
  fmpq_mat_init(dy_data->hankel_mat, d/2+1, d/2+1);
  for (i=0; i<=1; i++) fmpq_mat_init(dy_data->hankel_dets[i], d+1, 1);
  fmpq_set_si(fmpq_mat_entry(dy_data->hankel_dets[0], 0, 0), d, 1);
  fmpq_set_si(fmpq_mat_entry(dy_data->hankel_dets[1], 0, 0), 1, 1);
  
  dy_data->wlen = 3*d+7;
  dy_data->w = _fmpz_vec_init(dy_data->wlen);
  dy_data->w2len = 4;
  dy_data->w2 = _fmpq_vec_init(dy_data->w2len);

  return(dy_data);
}

/* Static memory deallocation. */
void ps_static_clear(ps_static_data_t *st_data) {
  if (st_data == NULL) return;
  int d = st_data->d;
  fmpz_clear(st_data->q);
  _fmpq_vec_clear(st_data->f, d+1);
  _fmpz_vec_clear(st_data->modlist, d+1);
  _fmpz_vec_clear(st_data->binom_mat, (d+1)*(d+1));
  _fmpz_vec_clear(st_data->sum_mats, (d+1)*(d+1));
  _fmpz_vec_clear(st_data->eval_pm2_mats, 2*(d+1)*(d+1));
  free(st_data);
}

/* Dynamic memory deallocation. */
void ps_dynamic_clear(ps_dynamic_data_t *dy_data) {
  if (dy_data == NULL) return;
  int i, d = dy_data->d;
  _fmpz_vec_clear(dy_data->pol, d+1);
  _fmpz_vec_clear(dy_data->sympol, 2*d+3);
  _fmpz_vec_clear(dy_data->upper, d+1);
  fmpq_mat_clear(dy_data->power_sums);
  fmpq_mat_clear(dy_data->hankel_mat);
  for (i=0; i<=1; i++) fmpq_mat_clear(dy_data->hankel_dets[i]);
  _fmpz_vec_clear(dy_data->w, dy_data->wlen);
  _fmpq_vec_clear(dy_data->w2, dy_data->w2len);
  free(dy_data);
}

/*****
  Low-level flow control
*****/

/* Increment the current moving counter and update stored data to match. 
   If step is NULL it is interpreted as 1. */

inline void step_forward(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int n, fmpz_t step) {
  int k = st_data->d - n;
  fmpz *pol = dy_data->pol;
  fmpq *f = st_data->f + n;
  fmpq *t2q = dy_data->w2;
  fmpq *t;
  fmpq_mat_struct *hankel_dets;

  if (step == NULL) {
    fmpq_set(t2q, f);
    fmpz_add(pol+n, pol+n, st_data->modlist+n);
  }
  else {
    fmpq_mul_fmpz(t2q, f, step);
    fmpz_addmul(pol+n, step, st_data->modlist+n);
  }
  t = fmpq_mat_entry(dy_data->power_sums, k, 0);
  fmpq_sub(t, t, t2q);
  hankel_dets = dy_data->hankel_dets[0];
  t = fmpq_mat_entry(hankel_dets, k, 0);
  if (k > 1) fmpq_submul(t, fmpq_mat_entry(hankel_dets, k-2, 0), t2q);
  else fmpq_sub(t, t, t2q);
  hankel_dets = dy_data->hankel_dets[1];
  t = fmpq_mat_entry(hankel_dets, k, 0);
  if (k > 1) fmpq_addmul(t, fmpq_mat_entry(hankel_dets, k-2, 0), t2q);
  else fmpq_add(t, t, t2q);
}

/* The following is the key subroutine: given some initial coefficients, compute
   a lower and upper bound for the next coefficient. Return 1 iff the resulting
   interval is nonempty.

   The value of dy_data->pol+n is assumed to be correct modulo st_data->modlist+n.
*/
int set_range_from_power_sums(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int n) {
  /* Static data */

  int i, j, r, s;
  int d = st_data->d;
  int k = d - n;
  int force_squarefree = st_data->force_squarefree;
  fmpz *modulus = st_data->modlist + n;
  fmpz *q = st_data->q;
  int q_is_1 = fmpz_is_one(q);
  fmpq *f = st_data->f + n;

  /* Dynamic data */

  fmpz *pol = dy_data->pol + n;
  fmpq_mat_struct *power_sums = dy_data->power_sums;
  fmpq_mat_struct *hankel_dets;

  /* Pointers into persistent working memory */

  fmpz *lower = dy_data->w;
  fmpz *upper = lower+1;
  fmpz *t0z = lower+2;
  fmpz *t1z = lower+3;
  fmpz *t2z = lower+4; // Affected by change_by_sign
  fmpz *tpol = lower+5; // Length d+1
  fmpz *f0 = lower+d+6; // Length d+1
  fmpz *f1 = lower+2*d+7; // Length d

  fmpq *t0q = dy_data->w2;
  fmpq *t1q = t0q+1;
  fmpq *t2q = t0q+2; // Affected by change_by_sign
  fmpq *t3q = t0q+3; // Affected by change_by_sign

  /* Unallocated pointers */

  fmpz *tz;
  fmpq *t, *t1, *t2;
  fmpq_mat_t hankel_mat;

  /* Adjust lower and upper bounds within set_range_from_power_sums.
     This overwrites t2z, t2q, and (if val2 != NULL) also t3q.
     The values in val1 and val2 need not be canonicalized.
     The pair (val1, val2) stands for val1 + val2*sqrt(q);
     passing NULL for val2 is a faster variant of passing 0.
     Given that this value is a monic linear function of the k-th power sum, then:

     -- passing r = 0 imposes the condition g >= 0 (or g > 0 if force_squarefree != 0);
     -- passing r = 1 imposes the condition g <= 0 (or g < 0 if force_squarefree != 0);
     -- passing update = 0 means we are setting the initial bounds;
     -- passing update = 1 means we are updating previously set bounds.
  */

  inline void change_by_sign(int update, int r, const fmpq_t val1, const fmpq_t val2) {
    fmpq_div_raw(t2q, val1, f);
    if (val2 == NULL) {
      if (r ^ force_squarefree) fmpq_ceil(t2z, t2q);
      else fmpq_floor(t2z, t2q);
    } else {
      fmpq_div_raw(t3q, val2, f);
      if (r ^ force_squarefree) fmpq_ceil_quad(t2z, t2q, t3q, q);
      else fmpq_floor_quad(t2z, t2q, t3q, q);
    }
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

  /* If k>d, no further coefficients to bound. */

  if (k>d) return(1);

  /* If modulus==0, reduce the interval to [0]. */

  i = fmpz_is_zero(modulus);
  if (i) {
    fmpz_zero(lower);
    fmpz_zero(upper);
  }

  /* Update power_sums[k] using the Girard-Newton formula. */

  #define POW(x) fmpq_mat_entry(power_sums, x, 0)
  tz = fmpq_numref(POW(0));
  fmpz_set_ui(tz, k); // Temporary change to apply Girard-Newton
  fmpq_mat_fmpz_vec_mul(t0q, pol, k, power_sums);
  fmpz_neg(t0z, pol+k);
  fmpq_div_fmpz(POW(k), t0q, t0z);
  fmpz_set_ui(tz, d); // Change back to the correct value, needed for Chebyshev criterion

  /* Chebyshev criterion: the k-th symmetrized power sum must lie in [-2*d*q^(k/2), 2*d*q^(k/2)]. */

  fmpq_mat_fmpz_vec_mul(t0q, st_data->sum_mats+(d+1)*k, k+1, power_sums);
  if (q_is_1 || k == 1) fmpz_set_ui(t0z, 2*d);
  else {
    fmpz_pow_ui(t1z, q, k/2);
    fmpz_mul_si(t0z, t1z, 2*d);
  }
  if (q_is_1 || k%2 == 0) {
    fmpq_sub_fmpz(t1q, t0q, t0z);
    fmpq_add_fmpz(t0q, t0q, t0z);
    t = NULL; t1 = t0q; t2 = t1q;
  } else {
    fmpq_set_fmpz(t1q, t0z);
    t = t1q; t1 = t0q; t2 = t0q;
  }
  change_by_sign(i, 0, t1, t); // i == fmpz_is_zero(modulus)
  if (t != NULL) fmpz_neg(fmpq_numref(t), fmpq_numref(t));
  change_by_sign(i, 1, t2, t);

  /* Descartes criterion: the evaluations of the n-th derivative of pol at -2*sqrt(q), 2*sqrt(q)
     have the correct signs. */
  _fmpz_vec_dot(t0z, st_data->eval_pm2_mats+(d+1)*(2*k), pol, k+1);
  _fmpz_vec_dot(t1z, st_data->eval_pm2_mats+(d+1)*(2*k+1), pol, k+1);
  fmpz_set(fmpq_denref(t0q), pol+k);
  fmpz_set(fmpq_denref(t1q), pol+k);
  if (q_is_1) {
    fmpz_add(fmpq_numref(t0q), t0z, t1z);
    fmpz_sub(fmpq_numref(t1q), t0z, t1z);
    t = NULL; t1 = t0q; t2 = t1q;
  } else {
    fmpz_set(fmpq_numref(t0q), t0z);
    fmpz_set(fmpq_numref(t1q), t1z);
    t = t1q; t1 = t0q; t2 = t0q;
  }
  change_by_sign(1, 1, t1, t);
  if (t != NULL) fmpz_neg(fmpq_numref(t), fmpq_numref(t));
  change_by_sign(1, 1-(k%2), t2, t);
  if ((s = fmpz_cmp(lower, upper)) > 0) return(0);

  /* Hausdorff criterion: the relevant Hankel matrices have nonnegative determinant. */

  for (r=0; r<=1; r++) {
    if (k%2 == 1 && !q_is_1) continue; // Hankel matrix is not defined over Q
    s = k/2 + !(r == 1 && k%2 == 0);
    if (r == 0 || k%2 == 0)
      fmpq_mat_window_init(hankel_mat, dy_data->hankel_mat, 0, 0, s, s);
    for (i=0; i<s; i++)
      for (j=0; j<s; j++) {
        t = fmpq_mat_entry(hankel_mat, i, j);
        if (i > 0 && j+1 < s) // This is a repeat of a previously computed entry
          fmpq_set(t, fmpq_mat_entry(hankel_mat, i-1, j+1));
        else if (r == 0 && k%2 == 0)
          fmpq_set(t, POW(i+j));
        else if (r == 1 && k%2 == 0) {
          fmpq_mul_ui(t, POW(i+j), 4);
  	  if (!q_is_1) fmpq_mul_fmpz(t, t, q);
	  fmpq_sub(t, t, POW(i+j+2));
	} else {
          fmpq_mul_ui(t, POW(i+j), 2);
	  if (r == 0) fmpq_add(t, t, POW(i+j+1));
	  else fmpq_sub(t, t, POW(i+j+1));
	}
      } // Final value of t will be used again
    hankel_dets = dy_data->hankel_dets[r];
    t1 = fmpq_mat_entry(hankel_dets, k, 0);
    fmpq_mat_det(t1, hankel_mat);
    if (r == 1 || k%2 == 0) // Otherwise reuse the window size
      fmpq_mat_window_clear(hankel_mat);
    if (k > 1 && fmpq_sgn(t2 = fmpq_mat_entry(hankel_dets, k-2, 0)) > 0)
      fmpq_div_raw(t0q, t1, t2);
    else fmpq_set(t0q, t); // t was set in the for loop
    if (r == 1) fmpz_neg(fmpq_numref(t0q), fmpq_numref(t0q));
    change_by_sign(1, r, t0q, NULL);
    if ((s = fmpz_cmp(lower, upper)) > 0) return(0);
  }

  /* Compute the divided n-th derivative of pol, answer in tpol. */

  tz = st_data->binom_mat + (d+2)*n;
  for (i=0; i<=k; i++) fmpz_mul(tpol+i, tz+i, pol+i);

  /* Rolle criterion: tpol has all roots real.
    Note: we do not call change_by_sign hereafter, so it is now safe to assign to t2z. */

  #define TEST_ROOTS(x) _fmpz_poly_all_real_roots(tpol, k+1, f0, f1, force_squarefree, modulus, x)
  if (!s) { /* Handle the case upper == lower directly. */
    if (!TEST_ROOTS(lower)) return(0);
  } else {
    /* Find a single value where the Rolle criterion holds. */
    fmpz_sub(t0z, upper, lower);
    fmpz_add_ui(t0z, t0z, 1);
    r = fmpz_flog_ui(t0z, 2); // r = floor(log_2 (upper-lower+1))
    s = 0;
    do {
      fmpz_one_2exp(t2z, r);
      if (!r) fmpz_set(t0z, lower);
      else {
        fmpz_add(t0z, lower, t2z);
        fmpz_sub_ui(t0z, t0z, 1);
      }
      do {
        if (!(s = TEST_ROOTS(t0z))) fmpz_addmul_ui(t0z, t2z, 2);
      } while (!s && fmpz_cmp(t0z, upper) <= 0);
      if (!s) r--;
    } while (!s && r >= 0);
    if (!s) return(0);

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
  }

  /* Set the new upper bound. */
  fmpz_mul(upper, upper, modulus);
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
void ps_dynamic_split(ps_dynamic_data_t *dy_data, ps_dynamic_data_t *dy_data2) {
  if ((dy_data == NULL) || (dy_data2 == NULL) || (dy_data->flag <= 0) || dy_data2->flag) return;

  int i, j, d = dy_data->d, n = dy_data->n, ascend = dy_data->ascend, k = n + ascend;

  for (i=d; i>k; i--)
    if (fmpz_cmp(dy_data->pol+i, dy_data->upper+i) < 0) {
      /* Copy the current state of the donor to the donee process. */
      dy_data2->n = n;
      dy_data2->ascend = ascend;
      _fmpz_vec_set(dy_data2->pol+k, dy_data->pol+k, d+1-k);
      _fmpz_vec_set(dy_data2->upper+k, dy_data->upper+k, d+1-k);
      fmpq_mat_set(dy_data2->power_sums, dy_data->power_sums);
      for (j=0; j<=1; j++) fmpq_mat_set(dy_data2->hankel_dets[j], dy_data->hankel_dets[j]);
      
      /* Restrict the donee process to the current branch. */
      fmpz_set(dy_data2->upper+i, dy_data2->pol+i);
      /* Remove the current branch from the donor process. */
      dy_data->ascend = i - n;
      /* Make the donee process available for further splitting. */
      dy_data2->flag = 1;
      break;
    }
}


/* Top-level flow control: allow one process to run for up to max_steps iterations,
   or until it finds a polynomial to be returned, whichever comes first.

   Return value sent back in dy_data->flag:
   1: in process
   2: found a solution (returned in dy_data->sympol)
   0: tree exhausted
   -1: maximum number of nodes reached

   It is *not* threadsafe to run next_pol and dynamic_split simultaneously on the same dy_data.
*/

void next_pol(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int max_steps) {
  if (dy_data==NULL || !dy_data->flag) return; // No work assigned to this process
  dy_data->flag = 0; // Prevent work-stealing while this process is running

  int d = st_data->d;
  int n = dy_data->n;
  if (n>d) return; // This process is exhausted.

  int ascend = dy_data->ascend;
  long node_limit = st_data->node_limit;
  long node_count = dy_data->node_count;
  fmpz *pol = dy_data->pol;
  fmpz *upper = dy_data->upper;

  int i, j, flag = 1, count_steps = 0;

  while (flag == 1 && count_steps <= max_steps) {
    count_steps += 1;
    if (ascend) { // Ascend the tree and step forward as needed.
      n += ascend;
      if (n > d) flag = 0; // This process is exhausted.
      else {
	ascend = (fmpz_cmp(pol+n, upper+n) >= 0);
	if (!ascend) step_forward(st_data, dy_data, n, NULL);
      }
    } else if (n < 0) { // Return a solution.
      fmpz *sympol = dy_data->sympol;
      fmpz *t;
      fmpz *q = st_data->q;
      int q_is_1 = fmpz_is_one(q);

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
    } else { // Compute children of the current node.
      n -= 1;
      ascend = !set_range_from_power_sums(st_data, dy_data, n);
      if (ascend) { // Found a terminal node
	node_count += 1;
	if (node_limit != -1 && node_count >= node_limit) flag = -1;
      }
    }
  }
  dy_data->ascend = ascend;
  dy_data->n = n;
  dy_data->node_count = node_count;
  dy_data->flag = flag;
}
