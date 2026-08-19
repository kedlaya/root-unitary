/*
  Low-level code to exhaust over trees of Weil polynomials.
  For examples of parallel execution, see the Cython, C, and Rust wrappers.

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

#include <stdlib.h>
#include "flint-utils.c"
#include "power_sums.h"

/* Check for OpenMP at runtime. */
int num_threads() {
  #if defined(_OPENMP)
  return omp_get_max_threads();
  #endif
  return(1);
}

/*****
  Memory allocation and deallocation
*****/

/* Static memory allocation and initialization. */
ps_static_data_t *ps_static_init(int d, const fmpz_t q, const fmpz_t lead, const fmpz *modlist,
                                 int num_constraints, const fmpz *constraints,
                                 long node_limit, int force_squarefree) {
  int i, j;
  fmpz *k0, *pol;

  ps_static_data_t *st_data = (ps_static_data_t *)malloc(sizeof(ps_static_data_t));
  st_data->d = d;
  fmpz_init_set(st_data->q, q);
  st_data->q_is_1 = fmpz_is_one(q);
  st_data->q_is_square = fmpz_is_square(q);
  st_data->lead_is_1 = fmpz_is_one(lead);
  if (st_data->q_is_square) {
     fmpz_init(st_data->q_sqrt);
     fmpz_sqrt(st_data->q_sqrt, q);
  }
  st_data->node_limit = node_limit;
  st_data->force_squarefree = force_squarefree;

  st_data->modlist = _fmpz_vec_init(d+1);
  k0 = st_data->modlist; // Used as a temporary variable for now
  
  fmpz_init(st_data->c0); // c0 = 4*lead^2*q
  fmpz_mul(st_data->c0, lead, lead);
  fmpz_mul(st_data->c0, st_data->c0, q);
  fmpz_mul_ui(st_data->c0, st_data->c0, 4);
  fmpz_init(st_data->c1); // c1 = 2*lead*sqrt(q)
  if (st_data->q_is_square) fmpz_sqrt(st_data->c1, st_data->c0);

  /* Matrix of cefficients of 2*(i-th Chebyshev polynomial)(x/2).
     The coefficient of x^j is multiplied by c^{i-j} q^{(i-j)/2} where c is the leading coefficient. */
  st_data->sum_mats = _fmpz_vec_init((d+1)*(d+1));
  for (i=0; i<=d; i++) {
    pol = st_data->sum_mats + (d+1)*i;
    _fmpz_poly_chebyshev_t(pol, i);
    _fmpz_poly_scale_2exp(pol, i+1, -1);
    if (!st_data->q_is_1)
      for (j=i%2; j<=i; j+=2) {
        fmpz_pow_ui(k0, q, (i-j)/2);
        fmpz_mul(pol+j, k0, pol+j);
      }
    if (!st_data->lead_is_1)
      for (j=i%2; j<=i; j+=2) {
        fmpz_pow_ui(k0, lead, i-j);
        fmpz_mul(pol+j, k0, pol+j);
      }
  }

  /* Moduli constraints. If not specified, fix the leading coefficient and impose no other constraints. */
  if (modlist == NULL) for (i=0; i<d; i++) fmpz_one(st_data->modlist+i);
  else for (i=0; i<=d; i++) fmpz_set(st_data->modlist+i, modlist+d-i);

  /* Matrix of binomial coefficients. */
  st_data->binom_mat = _fmpz_vec_init((d+1)*(d+1));
  for (i=0; i<=d; i++)
    for (j=0; j<=d; j++)
      fmpz_bin_uiui(st_data->binom_mat+(d+1)*i+j, j, i);

  /* Matrices to evaluate derivatives at +-2*sqrt(q). */
  st_data->eval_pm2_mats = _fmpz_vec_init(2*(d+1)*(d+1));
  for (i=0; i<=d; i++)
    for (j=0; j<=i; j++) {
      k0 = st_data->eval_pm2_mats + (d+1)*(2*i+j%2) + j;
      if (st_data->q_is_1) fmpz_one(k0);
      else if (st_data->q_is_square) fmpz_pow_ui(k0, st_data->q_sqrt, j);
      else fmpz_pow_ui(k0, q, j/2);
      fmpz_mul_2exp(k0, k0, j);
      fmpz_mul_si(k0, k0, -i);
      fmpz_mul(k0, k0, st_data->binom_mat + (d+2)*(d-i) + j);
    }

  st_data->ranges = _fmpz_vec_init(d+1);
  for (i=1; i<=d; i++) {
    k0 = st_data->ranges + i;
    if (st_data->q_is_square) fmpz_pow_ui(k0, st_data->q_sqrt, i);
    else fmpz_pow_ui(k0, q, i/2);
    fmpz_mul_ui(k0, k0, 2*d);
    for (j=0; j<i; j++) fmpz_mul(k0, k0, lead);
  }
  
  /* Matrix to compute reciprocal transform. */
  fmpz_mat_init(st_data->pol_to_sym, 2*d+1, d+1);
  fmpz_mat_zero(st_data->pol_to_sym);
  for (i=0; i<=d; i++) {
    k0 = fmpz_mat_entry(st_data->pol_to_sym, d-i, i);
    fmpz_one(k0);
    for (j=i; j>0; j--) {
      if (j==i) fmpz_set(fmpz_mat_entry(st_data->pol_to_sym, d+i, i), k0); 
      else fmpz_add(fmpz_mat_entry(st_data->pol_to_sym, d-i+2*j, i), fmpz_mat_entry(st_data->pol_to_sym, d-i+2*j, i), k0);
      if (!st_data->q_is_1) fmpz_mul(k0, k0, q);
      fmpz_mul_ui(k0, k0, j);
      fmpz_divexact_ui(k0, k0, i-j+1);
    }
  }
  
  if (!st_data->lead_is_1) {
    st_data->lead_pows = _fmpz_vec_init(d+1);
    for (i=0; i<=d; i++) fmpz_pow_ui(st_data->lead_pows+i, lead, i);
  }

  /* Linear conditions on asymmetrized coefficients. */
  st_data->num_constraints = num_constraints;
  if (num_constraints) {
    st_data->constraints = _fmpz_vec_init((d+1) * num_constraints);
    _fmpz_vec_set(st_data->constraints, constraints, (d+1)*num_constraints);
    st_data->constraint_lens = malloc(sizeof(int *) * num_constraints);
    for (i=0; i<num_constraints; i++) {
      for (j=d; j>0 && fmpz_is_zero(constraints+(d+1)*i+j); j--) {}
      if (j > 0) st_data->constraint_lens[i] = j;
    }
  }

  return(st_data);
}

/* Dynamic memory allocation and initialization.
   Call with coefflist == NULL to prepare an inactive process. */
ps_dynamic_data_t *ps_dynamic_init(int d, fmpz *coefflist) {
  ps_dynamic_data_t *dy_data = (ps_dynamic_data_t *)malloc(sizeof(ps_dynamic_data_t));

  dy_data->d = d;
  dy_data->n = d;
  dy_data->node_count = 0;
  dy_data->ascend = 0;
  dy_data->pol = _fmpz_vec_init(d+1);
  dy_data->sympol = _fmpz_vec_init(2*d+1);
  dy_data->upper = _fmpz_vec_init(d+1);
  dy_data->power_sums_num = _fmpz_vec_init(d+1);
  fmpz_set_si(dy_data->power_sums_num, d);
  _fmpz_vec_zero(dy_data->power_sums_num+1, d);

  dy_data->hankel_dets = _fmpz_vec_init(2*d+2);
  fmpz_set_si(dy_data->hankel_dets, d);
  fmpz_set_si(dy_data->hankel_dets+1, 1);

  dy_data->wlen = 3*d+9;
  dy_data->w = _fmpz_vec_init(dy_data->wlen);

  if (coefflist != NULL) {
    dy_data->flag = 1; // Activate this process
    _fmpz_vec_set(dy_data->pol, coefflist, d+1);
  } else dy_data->flag = 0; // Flag this process as inactive
  return(dy_data);
}

/* Static memory deallocation. */
void ps_static_clear(ps_static_data_t *st_data) {
  if (st_data == NULL) return;

  int d = st_data->d;

  if (fmpz_is_square(st_data->q)) fmpz_clear(st_data->q_sqrt);
  fmpz_clear(st_data->q);
  fmpz_clear(st_data->c0);
  fmpz_clear(st_data->c1);
  _fmpz_vec_clear(st_data->modlist, d+1);
  _fmpz_vec_clear(st_data->binom_mat, (d+1)*(d+1));
  _fmpz_vec_clear(st_data->sum_mats, (d+1)*(d+1));
  _fmpz_vec_clear(st_data->eval_pm2_mats, 2*(d+1)*(d+1));
  _fmpz_vec_clear(st_data->ranges, d+1);
  fmpz_mat_clear(st_data->pol_to_sym);
  if (!st_data->lead_is_1) _fmpz_vec_clear(st_data->lead_pows, d+1);
  if (st_data->num_constraints) {
    free(st_data->constraint_lens);
    _fmpz_vec_clear(st_data->constraints, (d+1)*st_data->num_constraints);
  }
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

void step_forward(const ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int n, fmpz_t step) {
  int d = st_data->d;
  int k = d - n;
  fmpz *pol = dy_data->pol;
  fmpz *poln = pol+n;
  fmpz *modulus = st_data->modlist + n;
  int modulus_is_1 = fmpz_is_one(modulus);
  fmpz *pow_num = dy_data->power_sums_num + k;
  fmpz *det = dy_data->hankel_dets + 2*k;
  fmpz *t0z = dy_data->w;

  if (!st_data->lead_is_1) {
    fmpz_pow_ui(t0z, pol+d, k-1);
    fmpz_mul_ui(t0z, t0z, k);
  } else fmpz_set_ui(t0z, k);
  fmpz_mul(t0z, t0z, modulus);
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

int ascend_step_forward(const ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int n) {
  int d = st_data->d;
  fmpz *pol = dy_data->pol;
  fmpz *upper = dy_data->upper;
  do {n++; } while ((n <= d) && (fmpz_cmp(pol+n, upper+n) >= 0));
  if (n <= d) step_forward(st_data, dy_data, n, NULL);
  return(n);
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
  fmpz *f0 = w+3; // Length k-1
  fmpz *f1 = w+k+2; // Length k

  #define TEST_ROOTS(x) _fmpz_poly_all_real_roots(pol, k, f0, f1, force_squarefree, x, modulus)

  /* Handle the case upper == lower directly. */
  fmpz_sub(t0z, upper, lower);
  if (fmpz_is_zero(t0z)) return(TEST_ROOTS(lower));

  /* Look for a single value where the Rolle criterion holds. */
  fmpz_add_ui(t0z, t0z, 1);
  r = fmpz_flog_ui(t0z, 2); // r = floor(log_2 (upper-lower+1)); forced to be positive
  fmpz_one_2exp(t2z, r);
  while (1) {
    if (r) {
      fmpz_add(t0z, lower, t2z);
      fmpz_sub_ui(t0z, t0z, 1);
    } else fmpz_set(t0z, lower);
    do {
      if ((s = TEST_ROOTS(t0z))) break; 
      else fmpz_addmul_ui(t0z, t2z, 2);
    } while (fmpz_cmp(t0z, upper) <= 0);
    if (s) break;
    if (--r < 0) return(0); // Found nothing
    fmpz_divexact_ui(t2z, t2z, 2);
  }

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
  while (!fmpz_equal(lower, t0z)) {
    fmpz_fmid(t2z, lower, t0z);
    if (TEST_ROOTS(t2z)) fmpz_set(t0z, t2z);
    else fmpz_add_ui(lower, t2z, 1);
  }
  while (!fmpz_equal(t1z, upper)) {
    fmpz_cmid(t0z, t1z, upper);
    if (TEST_ROOTS(t0z)) fmpz_set(t1z, t0z);
    else fmpz_sub_ui(upper, t0z, 1);
  }
  return(1);
}

/* Update pow_num[k] using the Girard-Newton formula.
   This is the k-th power sum times c^k where c is the leading coefficient. */

void update_power_sums(fmpz *pow_num, fmpz *pol, int k) {
  fmpz *tz = pow_num + k;
  fmpz *lead = pol + k;
  int i;

  if (fmpz_is_one(lead)) {
    tz = pow_num + k;
    _fmpz_vec_dot_general(tz, NULL, 1, pol+1, pow_num+1, 0, k-1);
    fmpz_submul_ui(tz, pol, k);
  } else {
    fmpz_mul_si(tz, pol, -k);
    for (i=1; i<k; i++) fmpz_fmms_wrapper(tz, tz, lead, pol+i, pow_num+i);
  }
}

/* The following is the key subroutine: given some initial coefficients, compute
   a lower and upper bound for the next coefficient. Return 1 iff the resulting
   interval is nonempty.

   The value of dy_data->pol+n is assumed to be correct modulo st_data->modlist+n.
*/
int set_range_from_power_sums(const ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int n) {
  if (n < 0) return(1);

  /* Static data */

  int d = st_data->d;
  int k = d - n;
  int k2 = k % 2;
  int force_squarefree = st_data->force_squarefree;
  fmpz *modulus = st_data->modlist + n;
  int modulus_is_0 = fmpz_is_zero(modulus);
  int modulus_is_1 = fmpz_is_one(modulus);
  const fmpz *q = st_data->q;
  int q_is_square = st_data->q_is_square;
  fmpz *lead_pow = (st_data->lead_is_1) ? NULL : st_data->lead_pows+k-1;

  /* Dynamic data */

  fmpz *pol = dy_data->pol + n;
  fmpz *pow_num = dy_data->power_sums_num;

  /* Integers allocated from working space, maintained throughout */

  fmpz *f = dy_data->w;
  fmpz *upper = f+1;
  fmpz *lower = f+2;

  /* Local working space, not maintained throughout */

  int i, j, j1, s;
  fmpz *w = f+3; // Length 3*k+3, passed to subroutines

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
    fmpq_floor_ceil_quad(t2z, r ^ force_squarefree, val1_num, val1_den, val2_num, f, q);
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

  update_power_sums(pow_num, pol, k);

  /* Chebyshev criterion: the k-th symmetrized power sum must lie in [-2*d*q^(k/2), 2*d*q^(k/2)]. */

  _fmpz_vec_dot(t0z, st_data->sum_mats+(d+1)*k, pow_num, k+1);
  fmpz_set(t1z, st_data->ranges+k);
  if (q_is_square || k2 == 0) { // q^{k/2} is rational
    fmpz_add(t2z, t0z, st_data->ranges+k);
    change_by_sign(modulus_is_0, 0, t2z, lead_pow, NULL);
    fmpz_sub(t2z, t0z, st_data->ranges+k);
    change_by_sign(modulus_is_0, 1, t2z, lead_pow, NULL);
  } else { // q^{k/2} is irrational
    change_by_sign(modulus_is_0, 0, t0z, lead_pow, st_data->ranges+k);
    fmpz_neg(t1z, st_data->ranges+k);
    change_by_sign(modulus_is_0, 1, t0z, lead_pow, t1z);
  }

  /* Impose optional linear constraints on power sums. */

  for (i=0; i<st_data->num_constraints; i++)
    if (st_data->constraint_lens[i] == k) {
      tz = st_data->constraints+(d+1)*i;
      _fmpz_vec_dot(t0z, tz, pow_num, k+1);
      tz += k;
      if (fmpz_sgn(tz) < 0) {
        fmpz_neg(t0z, t0z);
        fmpz_neg(t1z, tz);
        change_by_sign(1, 1, t0z, t1z, NULL);
      } else change_by_sign(1, 0, t0z, tz, NULL);
    }

  /* Descartes criterion: the evaluations of the n-th derivative of pol at -2*sqrt(q), 2*sqrt(q)
     have the correct signs. */

  _fmpz_vec_dot(t0z, st_data->eval_pm2_mats+(d+1)*(2*k), pol, k+1);
  _fmpz_vec_dot(t1z, st_data->eval_pm2_mats+(d+1)*(2*k+1), pol, k+1);
  if (q_is_square) {
    fmpz_add(t2z, t0z, t1z);
    change_by_sign(1, 1, t2z, NULL, NULL);
    fmpz_sub(t2z, t0z, t1z);
    change_by_sign(1, 1-k2, t2z, NULL, NULL);
  } else {
    change_by_sign(1, 1, t0z, NULL, t1z);
    fmpz_neg(t1z, t1z);
    change_by_sign(1, 1-k2, t0z, NULL, t1z);
  }
  if (fmpz_cmp(lower, upper) > 0) return(0);

  /* Hausdorff criterion: the relevant Hankel matrices have nonnegative determinant.
     In order to restrict to integer arithmetic, we skip this condition when
     k is odd and q is not a perfect square. */

  if (q_is_square || k2 == 0) // Hankel matrix is defined over Q
    for (i=1; i>=0; i--) {
      s = (k2 == 1) ? k : (i == 0) ? k + 1 : k - 1;
      tz = dy_data->hankel_dets + 2*k + i;
      j = k > 1 && (force_squarefree || fmpz_sgn(tz-4) > 0);
      if (q_is_square && i == 0 && k > 1 && (force_squarefree || fmpz_sgn(tz-3) > 0)) {
        /* Use the recursive relationship between upper and lower Hankel determinants. */
        if (k2 == 1) {
           tza = t1z;
           fmpz_mul_ui(tza, tz-2, 2);
           fmpz_mul(tza, tza, st_data->c1);
        } else tza = tz-2;
        fmpz_fmms_wrapper(t1z, tz-1, tza, tz-4, tz+1);
        fmpz_divexact(tz, t1z, tz-3);
        if (k2 == 1 && !j) { // Need the corner entry for a linear constraint
          tza = w;
          fmpz_mul(tza+s-1, pow_num+s-1, st_data->c1);
          fmpz_add(tza+s-1, tza+s-1, pow_num+s+1-k2);
        } else tza = pow_num;
      } else {
        /* Compute the Hankel determinant by condensation (if possible) or directly. 
           Predefined: c1 == 2*lead*sqrt(q), c0 == 4*lead^2*q. */
        if (i == 1 || k2 == 1) {
          tza = (fmpz *)(k2 == 1 ? st_data->c1 : st_data->c0);
          j1 = k2 == 0 ? 2 : 1;
          if (i == 0) _fmpz_vec_scalar_fmma_one(w, pow_num, tza, pow_num+j1, s);
          else _fmpz_vec_scalar_fmms_one(w, pow_num, tza, pow_num+j1, s);
          tza = w;
        } else tza = pow_num;
        if (!hankel_determinant_condensation(tz, tza, s, w+k+1))
          hankel_determinant_direct(tz, tza, s);
      }

      /* Deduce a linear constraint. */
      if (j) {
        if (lead_pow != NULL) { tza = t1z; fmpz_mul(tza, tz-4, lead_pow); }
        else tza = tz-4;
        if (i == 1) { fmpz_neg(t0z, tz); tz = t0z; }
      } else { // Argue that the corner entry is nonnegative
        if (i == 1) { tz = t0z; fmpz_neg(tz, tza+s-1); }
        else tz = tza+s-1;
        tza = lead_pow;
      }
      change_by_sign(1, i, tz, tza, NULL);
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
   Return value is 1 if a split is executed and 0 otherwise.

   It is safe to run this in parallel as long as the instances of dy_data are pairwise distinct
   and the instances of dy_data2 are pairwise distinct.
*/
int ps_dynamic_split(const ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, ps_dynamic_data_t *dy_data2) {
  if (dy_data == NULL || dy_data2 == NULL || dy_data == dy_data2 || dy_data->flag <= 0 || dy_data2->flag) return(0);

  int i;
  int d = dy_data->d;
  int k = dy_data->n + dy_data->ascend;

  for (i=d; i>k; i--)
    if (fmpz_cmp(dy_data->pol+i, dy_data->upper+i) < 0) {
      int j = d - i + 1;
      fmpz *lower = dy_data->pol+i;
      fmpz *upper = dy_data->upper+i;
      fmpz *modulus = st_data->modlist+i;
      int modulus_is_1 = fmpz_is_one(modulus);
      fmpz *t0z = dy_data->w; // Temporary variable

      /* Copy the current state of the donor to the donee process. */
      _fmpz_vec_set(dy_data2->pol+i, lower, j);
      _fmpz_vec_set(dy_data2->upper+i, upper, j);
      _fmpz_vec_set(dy_data2->power_sums_num, dy_data->power_sums_num, j);
      _fmpz_vec_set(dy_data2->hankel_dets, dy_data->hankel_dets, 2*j);

      /* Clear unused values to make sure large mpz's get deallocated promptly. */
      _fmpz_vec_zero(dy_data2->power_sums_num+j, d+1-j);
      _fmpz_vec_zero(dy_data2->hankel_dets+2*j, 2*(d+1-j));
      _fmpz_vec_zero(dy_data->w, dy_data->wlen);
      _fmpz_vec_zero(dy_data2->w, dy_data2->wlen);

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
      return(1);
    }
  return(0);
}

/* Compute the reciprocal transform of pol and store it in sympol. */

void reciprocal_transform(const ps_static_data_t *st_data, ps_dynamic_data_t *dy_data) {
  int d = st_data->d;
  fmpz *pol = dy_data->pol;
  fmpz *sympol = dy_data->sympol;

  fmpz_mat_mul_fmpz_vec(sympol, st_data->pol_to_sym, pol, d+1);
}

/* Top-level flow control: allow one process to run for up to max_steps iterations,
   or until it finds a polynomial to be returned, whichever comes first.

   Return value also sent back in dy_data->flag:
   1: in process
   2: found a solution (returned in dy_data->sympol)
   0: tree exhausted
   -1: maximum number of nodes reached

   It is threadsafe to run this in parallel with dynamic_split.
*/

int ps_next_pol(const ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int max_steps) {
  if (dy_data==NULL || !dy_data->flag) return(0); // No work assigned to this process

  int d = st_data->d;
  int n = dy_data->n;
  if (n>d) return(0); // This process is exhausted.

  int ascend = dy_data->ascend;
  long node_limit = st_data->node_limit;
  long node_count = dy_data->node_count;

  int flag = 1;
  int count_steps = 0;
  int i;

  while (1) {
    if (ascend) { // Ascend the tree and step forward as needed.
      n = ascend_step_forward(st_data, dy_data, n);
      if (n > d) {flag = 0; break;}
      ascend = 0;
    } else if (n < 0) { // Return a solution.
      reciprocal_transform(st_data, dy_data);
      ascend = 1;
      flag = 2;
      break;
    } else { // Compute children of the current node.
      n--;
      if ((ascend = !set_range_from_power_sums(st_data, dy_data, n))) // Found a terminal node
	if (node_limit != -1 && (++node_count) >= node_limit) {
	  flag = -1;
	  break;
	}
      i = d-n+1; count_steps += i*i;
      if (count_steps > max_steps) break;
    }
  }
  /* Record the final working state. */
  dy_data->ascend = ascend;
  dy_data->n = n;
  dy_data->node_count = node_count;
  dy_data->flag = flag;
  return(flag);
}

void ps_cleanup(int n) {
  if (n == 0) flint_cleanup_master();
  else if (n == 1) flint_cleanup();
}
