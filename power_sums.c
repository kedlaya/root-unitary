/*
  Low-level code to exhaust over trees of Weil polynomials.
  This code does not implement parallelism; see the Cython wrapper.

  TODO: check for memory leaks.
  TODO: try the Routh-Hurwitz criterion.

#*****************************************************************************
#       Copyright (C) 2019 Kiran S. Kedlaya <kskedl@gmail.com>
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

/*****
  Arithmetic functions
  
  Note: input fmpq's need *not* be canonicalized, and output is not guaranteed to be either.
*****/

/* Set res to -a without canonicalizing. */
inline void fmpq_neg_raw(fmpq_t res, fmpq_t a) {
  fmpz_neg(fmpq_numref(res), fmpq_numref(a));
  fmpz_set(fmpq_denref(res), fmpq_denref(a));
}

/* Set res to a-b without canonicalizing. */
inline void fmpq_sub_raw(fmpq_t res, fmpq_t a, fmpq_t b) {
 fmpz *anum = fmpq_numref(a);
 fmpz *aden = fmpq_denref(a);
 fmpz *bnum = fmpq_numref(b);
 fmpz *bden = fmpq_denref(b);
 fmpz_fmms(fmpq_numref(res), anum, bden, bnum, aden);
 fmpz_mul(fmpq_denref(res), aden, bden);
}

/* Set res to a*b without canonicalizing. */
inline void fmpq_mul_raw(fmpq_t res, fmpq_t a, fmpq_t b) {
  fmpz_mul(fmpq_numref(res), fmpq_numref(a), fmpq_numref(b));
  fmpz_mul(fmpq_denref(res), fmpq_denref(a), fmpq_denref(b));
}

/* Set res to a/b without canonicalizing. */
inline void fmpq_div_raw(fmpq_t res, fmpq_t a, fmpq_t b) {
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
  fmpz_fdiv_q_2exp(res, res, 1);
}

/* Set res to ceil((a+b)/2). */
inline void fmpz_cmid(fmpz_t res, const fmpz_t a, const fmpz_t b) {
  fmpz_add(res, a, b);
  fmpz_cdiv_q_2exp(res, res, 1);
}

/* Set res to floor(sqrt(a)). */
inline void fmpz_sqrt_f(fmpz_t res, const fmpz_t a) {
  fmpz_sqrt(res, a);
}

/* Set res to ceil(sqrt(a)). */
inline void fmpz_sqrt_c(fmpz_t res, const fmpz_t a) {
  fmpz_sqrt(res, a);
  if (!fmpz_is_square(a)) fmpz_add_ui(res, res, 1);
}

/* Set res to floor(a + b sqrt(q)). */
inline void fmpq_floor_quad(fmpz_t res, fmpq_t a,
		     fmpq_t b, const fmpz_t q, int q_is_1) {
  if (b==NULL) fmpq_floor(res, a);
  else {
    fmpz *anum = fmpq_numref(a);
    fmpz *aden = fmpq_denref(a);
    fmpz *bnum = fmpq_numref(b);
    fmpz *bden = fmpq_denref(b);
    int bden_s = fmpz_sgn(bden);

    if (q_is_1) 
      fmpz_fmma(res, anum, bden, bnum, aden);
    else {
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
    }
    if (bden_s > 0) fmpz_fdiv_q(res, res, aden);
    else fmpz_cdiv_q(res, res, aden);
    fmpz_fdiv_q(res, res, bden);
  }
}

/* Set res to ceil(a + b sqrt(q)). */
inline void fmpq_ceil_quad(fmpz_t res, fmpq_t a,
		     fmpq_t b, const fmpz_t q, int q_is_1) {
  if (b==NULL) fmpq_ceil(res, a);
  else {
    fmpz *anum = fmpq_numref(a);
    fmpz *aden = fmpq_denref(a);
    fmpz *bnum = fmpq_numref(b);
    fmpz *bden = fmpq_denref(b);
    int bden_s = fmpz_sgn(bden);

    if (q_is_1) 
      fmpz_fmma(res, anum, bden, bnum, aden);
    else {
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
    }
    if (bden_s > 0) fmpz_cdiv_q(res, res, aden);
    else fmpz_fdiv_q(res, res, aden);
    fmpz_cdiv_q(res, res, bden);
  }
}

/*
    Use a subresultant (Sturm-Habicht) sequence to test whether a given
    polynomial has all real roots. Note that this test has an early abort
    mechanism: having all real roots means that the sign sequence has
    the maximal number of sign changes, so the test aborts as soon
    as a sign change is missed.

    This function assumes that:
        - {poly, n} is a normalized vector with n >= 2
        - {w, 2*n+1} is scratch space.
    If a and b are not NULL, we add a*b to the constant term before testing.

    Based on code by Sebastian Pancratz from the FLINT repository.
    TODO: reimplement using a half-GCD algorithm.
*/

int _fmpz_poly_all_real_roots(fmpz *poly, long n, fmpz *w, int force_squarefree,
			      const fmpz_t a, const fmpz_t b) {
  fmpz *f0     = w + 0*n;
  fmpz *f1     = w + 1*n;
  fmpz *c      = w + 2*n;
  fmpz *d      = w + 2*n+1;
  fmpz *t; // Not allocated, only used to swap pointers

  _fmpz_vec_set(f0, poly, n);
  /* Sanitize input so that n = deg(f0). */
  while ((n > 2) && fmpz_is_zero(f0+n-1))
    n--;
  if (n <= 2) return(1);
  if (a != NULL && b != NULL) fmpz_addmul(f0, a, b);
  _fmpz_poly_derivative(f1, f0, n);
  n--;
  int sgn0_l = fmpz_sgn(f0+n);

  while (1) {
    /* At this point deg(f0) = n, deg(f1) = n-1.
       We compute the negated pseudoremainder of f0 modulo f1 in two steps:
       f0 --> f1[n-1]*f0 - f0[n]*x*f1
       f0 --> f0[n-1]*f1 - f1[n-1]*f0

       f0 --> f1[n-1]*f0 - f0[n]*x*f1
       f0 --> f0[n-1]*f1 - f1[n-1]*f0
    */
    fmpz_set(c, f0+n);
    _fmpz_vec_scalar_mul_fmpz(f0, f0, n, f1+n-1);
    _fmpz_vec_scalar_submul_fmpz(f0+1, f1, n-1, c);
    n--;
    fmpz_set(c, f0+n);
    fmpz_neg(d, f1+n);
    _fmpz_vec_scalar_mul_fmpz(f0, f0, n, d);
    _fmpz_vec_scalar_addmul_fmpz(f0, f1, n, c);

    /* If f0 = 0, we win unless we are insisting on squarefree. */
    if (!force_squarefree && _fmpz_vec_is_zero(f0, n)) return(1);

    /* If we miss any one sign change, we cannot have enough. */
    if (fmpz_sgn(f0+n-1) != sgn0_l) return(0);

    /* If f0 is a scalar, it is nonzero and we win. */
    if (n == 1) return(1);

    /* Extract content from f0.
       This seems to do better in practice than an explicit subresultant computation. */
    _fmpz_vec_content(c, f0, n);
    _fmpz_vec_scalar_divexact_fmpz(f0, f0, n, c);

    /* Swap f0 with f1 at the pointer level. */
    t = f0; f0 = f1; f1 = t;
  }
}

/*****
  High-level memory management
*****/

/* Memory allocation and initialization. */
ps_static_data_t *ps_static_init(int d, fmpz_t q, int coeffsign, fmpz_t lead,
				 fmpz *modlist, long node_limit, int force_squarefree) {
  int i, j, k, l;
  ps_static_data_t *st_data;
  fmpz_poly_t pol;
  fmpz_t m;
  fmpz *k0;

  fmpz_poly_init(pol);
  fmpz_init(m);

  st_data = (ps_static_data_t *)malloc(sizeof(ps_static_data_t));

  st_data->d = d;
  st_data->sign = coeffsign;
  fmpz_init_set(st_data->q, q);
  st_data->node_limit = node_limit;
  st_data->force_squarefree = force_squarefree;

  fmpz_init_set(st_data->lead, lead);

  st_data->modlist = _fmpz_vec_init(d+1);
  st_data->f = _fmpq_vec_init(d+1);
  for (i=0; i<=d; i++) {
    fmpz_set(st_data->modlist+i, modlist+d-i);
    fmpq_set_si(st_data->f+i, d-i, 1);
    fmpq_div_fmpz(st_data->f+i, st_data->f+i, st_data->lead);
    /* In order to apply power sums and Descartes' rule of signs
       when the modulus is 0, we must pretend that the modulus is 1. */
    if (!fmpz_is_zero(st_data->modlist+i))
      fmpq_mul_fmpz(st_data->f+i, st_data->f+i, st_data->modlist+i);
  }

  fmpz_mat_init(st_data->binom_mat, d+1, d+1);
  for (i=0; i<=d; i++)
    for (j=0; j<=d; j++)
      fmpz_bin_uiui(fmpz_mat_entry(st_data->binom_mat, i, j), i, j);

  st_data->sum_mats = _fmpz_vec_init((d+1)*(d+1));
  _fmpz_vec_zero(st_data->sum_mats, (d+1)*(d+1)); // Redundant?
  for (i=0; i<=d; i++) {
    arith_chebyshev_t_polynomial(pol, i);
    for (j=i%2; j<=i; j+=2) {
      /* Coefficients of 2*(i-th Chebyshev polynomial)(x/2).
         If q != 1, the coeff of x^j is multiplied by q^{floor(i-j)/2}. */
      k0 = st_data->sum_mats+(d+1)*i+j;
      fmpz_set(k0, fmpz_poly_get_coeff_ptr(pol, j));
      if (j == 0) fmpz_mul_2exp(k0, k0, 1);
      else fmpz_fdiv_q_2exp(k0, k0, j-1);
      if (!fmpz_is_one(st_data->q)) {
        fmpz_pow_ui(m, st_data->q, (i-j)/2);
        fmpz_mul(k0, k0, m);
      }
    }
  }

  fmpz_poly_clear(pol);
  fmpz_clear(m);

  st_data->eval_pm2_mats = _fmpz_vec_init(2*(d+1)*(d+1));
  _fmpz_vec_zero(st_data->eval_pm2_mats, 2*(d+1)*(d+1)); // Redundant?
  for (i=0; i<=d; i++) {
    for (j=0; j<=i; j++) {
      k0 = st_data->eval_pm2_mats+(d+1)*(2*i+j%2)+j;
      fmpz_pow_ui(k0, st_data->q, j/2);
      fmpz_mul_2exp(k0, k0, j);
      fmpz_mul_si(k0, k0, -i);
    }
  }

  return(st_data);
}

ps_dynamic_data_t *ps_dynamic_init(int d, fmpz_t q, fmpz *coefflist) {
  ps_dynamic_data_t *dy_data;
  int i;

  dy_data = (ps_dynamic_data_t *)malloc(sizeof(ps_dynamic_data_t));
  dy_data->d = d;
  dy_data->q_is_1 = fmpz_is_one(q);

  /* Initialize mutable quantities */
  dy_data->n = d;
  dy_data->node_count = 0;
  dy_data->ascend = 0;
  dy_data->pol = _fmpz_vec_init(d+1);
  dy_data->sympol = _fmpz_vec_init(2*d+3);
  if (coefflist != NULL) {
    dy_data->flag = 1; // Activate this process
    for (i=0; i<=d; i++)
      fmpz_set(dy_data->pol+i, coefflist+i);
  } else dy_data->flag = 0;
  dy_data->upper = _fmpz_vec_init(d+1);

  fmpq_mat_init(dy_data->power_sums, d+1, 1);
  fmpq_set_si(fmpq_mat_entry(dy_data->power_sums, 0, 0), d, 1);
  fmpq_mat_init(dy_data->sum_prod, 1, 1);
  fmpq_mat_init(dy_data->hankel_mat, d/2+1, d/2+1);
  for (i=0; i<=1; i++) fmpq_mat_init(dy_data->hankel_dets[i], d+1, 1);
  fmpq_set_si(fmpq_mat_entry(dy_data->hankel_dets[0], 0, 0), d, 1);
  fmpq_set_si(fmpq_mat_entry(dy_data->hankel_dets[1], 0, 0), 1, 1);
  
  /* Allocate scratch space */
  dy_data->wlen = 3*d+10;
  dy_data->w = _fmpz_vec_init(dy_data->wlen);
  dy_data->w2len = 5;
  dy_data->w2 = _fmpq_vec_init(dy_data->w2len);
  return(dy_data);
}

/* Split off a subtree (without allocating new memory).
   The first process gives up on the current branch, up to the first coefficient that is not uniquely specified;
   the remaining work is yielded to the second process, which may in turn be split immediately.
*/
void ps_dynamic_split(ps_dynamic_data_t *dy_data, ps_dynamic_data_t *dy_data2) {
  if ((dy_data == NULL) || (dy_data->flag <= 0) || dy_data2->flag) return;

  int i, j, d = dy_data->d, n = dy_data->n, ascend = dy_data->ascend;

  for (i=d; i>n+ascend; i--)
    if (fmpz_cmp(dy_data->pol+i, dy_data->upper+i) <0) {
      dy_data2->n = n;
      dy_data2->ascend = ascend;
      _fmpz_vec_set(dy_data2->pol, dy_data->pol, d+1);
      _fmpz_vec_set(dy_data2->upper, dy_data->upper, d+1);
      fmpq_mat_set(dy_data2->power_sums, dy_data->power_sums);
      for (j=0; j<=1; j++) fmpq_mat_set(dy_data2->hankel_dets[j], dy_data->hankel_dets[j]);
      fmpz_set(dy_data2->upper+i, dy_data2->pol+i);
      dy_data->ascend = i-n;
      dy_data2->flag = 1; // This process can now itself be split.
      return;
  }
  return;
}

/* Memory deallocation. */
void ps_static_clear(ps_static_data_t *st_data) {
  if (st_data == NULL) return;
  int i, d = st_data->d;
  fmpz_clear(st_data->lead);
  fmpz_clear(st_data->q);
  fmpz_mat_clear(st_data->binom_mat);
  _fmpq_vec_clear(st_data->f, d+1);
  _fmpz_vec_clear(st_data->modlist, d+1);
  _fmpz_vec_clear(st_data->sum_mats, (d+1)*(d+1));
  _fmpz_vec_clear(st_data->eval_pm2_mats, 2*(d+1)*(d+1));
  free(st_data);
}

void ps_dynamic_clear(ps_dynamic_data_t *dy_data) {
  int i;
  
  if (dy_data == NULL) return;
  int d = dy_data->d;
  _fmpz_vec_clear(dy_data->pol, d+1);
  _fmpz_vec_clear(dy_data->sympol, 2*d+3);
  _fmpz_vec_clear(dy_data->upper, d+1);
  fmpq_mat_clear(dy_data->power_sums);
  fmpq_mat_clear(dy_data->sum_prod);
  fmpq_mat_clear(dy_data->hankel_mat);
  for (i=0; i<=1; i++) fmpq_mat_clear(dy_data->hankel_dets[i]);
  _fmpz_vec_clear(dy_data->w, dy_data->wlen);
  _fmpq_vec_clear(dy_data->w2, dy_data->w2len);
  free(dy_data);
}

/* Subroutines to adjust lower and upper bounds within set_range_from_power_sums.
   These use t0z, t0q, t4q as persistent scratch space.
   The pair (val1, val2) stands for val1 + val2*sqrt(q);
   passing NULL for val2 is a faster variant of passing 0.

  Usage: if g is a monic linear function of the k-th power sum, then
  set_upper(g) or change_upper(g) imposes the condition g >= 0;
  set_lower(g) or change_lower(g) imposes the condition g <= 0.

*/

#define STATE lower, upper, q, f, t0z, t0q, t4q, q_is_1
#define STATE_DECLARE fmpz_t lower, fmpz_t upper, fmpz_t q, fmpq_t f, fmpz_t t0z, fmpq_t t0q, fmpq_t t4q, int q_is_1
#define POL tpol, k+1, tpol2, st_data->force_squarefree, modulus

inline void set_lower(const fmpq_t val1, const fmpq_t val2, STATE_DECLARE) {
  fmpq_div_raw(t0q, val1, f);
  if (val2==NULL) fmpq_ceil(lower, t0q);
  else {
    fmpq_div_raw(t4q, val2, f);
    fmpq_ceil_quad(lower, t0q, t4q, q, q_is_1);
  }
}

inline void set_upper(const fmpq_t val1, const fmpq_t val2, STATE_DECLARE) {
  fmpq_div_raw(t0q, val1, f);
  if (val2==NULL) fmpq_floor(upper, t0q);
  else {
    fmpq_div_raw(t4q, val2, f);
    fmpq_floor_quad(upper, t0q, t4q, q, q_is_1);
  }
}

inline void change_lower(const fmpq_t val1, const fmpq_t val2, STATE_DECLARE) {
  fmpq_div_raw(t0q, val1, f);
  if (val2==NULL) fmpq_ceil(t0z, t0q);
  else {
    fmpq_div_raw(t4q, val2, f);
    fmpq_ceil_quad(t0z, t0q, t4q, q, q_is_1);
  }
  if (fmpz_cmp(t0z, lower) > 0) fmpz_set(lower, t0z);
}

inline void change_upper(const fmpq_t val1, const fmpq_t val2, STATE_DECLARE) {
  fmpq_div_raw(t0q, val1, f);
  if (val2==NULL) fmpq_floor(t0z, t0q);
  else {
    fmpq_div_raw(t4q, val2, f);
    fmpq_floor_quad(t0z, t0q, t4q, q, q_is_1);
  }
  if (fmpz_cmp(t0z, upper) < 0) fmpz_set(upper, t0z);
}

inline void change_by_sign(const fmpq_t val1, const fmpq_t val2, int s, STATE_DECLARE) {
  if (s == 0) change_upper(val1, val2, STATE);
  else {
    fmpq_neg_raw(t0q, val1);
    change_lower(t0q, val2, STATE);
  }
}

inline void change_lower_strict(const fmpq_t val1, const fmpq_t val2, STATE_DECLARE) {
  fmpq_div_raw(t0q, val1, f);
  if (val2==NULL) fmpq_floor(t0z, t0q);
  else {
    fmpq_div_raw(t4q, val2, f);
    fmpq_floor_quad(t0z, t0q, t4q, q, q_is_1);
  }
  fmpz_add_ui(t0z, t0z, 1);
  if (fmpz_cmp(t0z, lower) > 0) fmpz_set(lower, t0z);
}

inline void change_upper_strict(const fmpq_t val1, const fmpq_t val2, STATE_DECLARE) {
  fmpq_div_raw(t0q, val1, f);
  if (val2==NULL) fmpq_ceil(t0z, t0q);
  else {
    fmpq_div_raw(t4q, val2, f);
    fmpq_ceil_quad(t0z, t0q, t4q, q, q_is_1);
  }
  fmpz_sub_ui(t0z, t0z, 1);
  if (fmpz_cmp(t0z, upper) < 0) fmpz_set(upper, t0z);
}

/* The following is the key subroutine: given some initial coefficients, compute
   a lower and upper bound for the next coefficient. Return 1 iff the resulting
   interval is nonempty.

*/
int set_range_from_power_sums(ps_static_data_t *st_data,
			      ps_dynamic_data_t *dy_data) {
  int i, j, r, s;
  int d = st_data->d;
  int n = dy_data->n;
  int k = d+1-n;
  int q_is_1 = dy_data->q_is_1;
  fmpz *modulus = st_data->modlist+n-1;
  fmpz *pol = dy_data->pol;
  fmpz *q = st_data->q;
  fmpq *f = (fmpq *)(st_data->f+n-1);
  fmpq *t, *t1, *t2; // Unallocated, will be assigned from existing pointers
  fmpq_mat_t hausdorff_prod, hankel_mat;

  /* Allocate temporary variables from persistent scratch space. */
  fmpz *tpol = dy_data->w;
  fmpz *tpol2 = dy_data->w+d+1;

  fmpz *t0z = dy_data->w+3*d+5; // This gets overwritten by subroutines
  fmpz *t1z = dy_data->w+3*d+6;
  fmpz *t2z = dy_data->w+3*d+7;
  fmpz *lower = dy_data->w+3*d+8;
  fmpz *upper = dy_data->w+3*d+9;

  fmpq *t0q = dy_data->w2; // This gets overwritten by subroutines
  fmpq *t1q = dy_data->w2+1;
  fmpq *t2q = dy_data->w2+2;
  fmpq *t3q = dy_data->w2+3;
  fmpq *t4q = dy_data->w2+4; // This gets overwritten by subroutines

  /* If k>d, no further coefficients to bound. */
  if (k>d) return(1);

  /* Update power_sums[k] using the Girard-Newton formula. */
  t = fmpq_mat_entry(dy_data->power_sums, 0, 0);
  fmpq_set_si(t, k, 1);  
  fmpq_mat_fmpz_vec_mul(t0q, pol+d-k, k, dy_data->power_sums);
  fmpz_neg(t0z, pol+d);
  fmpq_div_fmpz(fmpq_mat_entry(dy_data->power_sums, k, 0), t0q, t0z);
  fmpq_set_si(t, d, 1);

  /* Condition: the k-th symmetrized power sum must lie in [-2*sqrt(q), 2*sqrt(q)]. */
  fmpq_mat_fmpz_vec_mul(t2q, st_data->sum_mats+(d+1)*k, k+1, dy_data->power_sums);
  fmpz_set_si(t1z, 2*d);
  if (!q_is_1) {
    fmpz_pow_ui(t0z, q, k/2);
    fmpz_mul(t1z, t1z, t0z);
  }
  if (k%2==0 || q_is_1) {
    fmpq_sub_fmpz(t0q, t2q, t1z);
    set_lower(t0q, NULL, STATE);
    fmpq_add_fmpz(t0q, t2q, t1z);
    set_upper(t0q, NULL, STATE);
  } else {
    fmpq_one(t1q);
    fmpz_set(fmpq_numref(t1q), t1z);
    set_upper(t2q, t1q, STATE);
    fmpq_neg(t1q, t1q);
    set_lower(t2q, t1q, STATE);
  }

  /* Compute the divided (n-1)-st derivative of pol, answer in tpol. */
  for (i=0; i<=k; i++)
    fmpz_mul(tpol+i, fmpz_mat_entry(st_data->binom_mat, n-1+i, n-1), pol+n-1+i);

  /* Condition: Descartes' rule of signs applies at -2*sqrt(q), +2*sqrt(q).
   This is only a new condition for the evaluations at these points. */

  _fmpz_vec_dot(fmpq_numref(t1q), st_data->eval_pm2_mats+(d+1)*(2*k), tpol, k+1);
  fmpz_set(fmpq_denref(t1q), pol+d);
  _fmpz_vec_dot(fmpq_numref(t2q), st_data->eval_pm2_mats+(d+1)*(2*k+1), tpol, k+1);
  fmpz_set(fmpq_denref(t2q), pol+d);

  /* If checking for squarefree, shear endpoints off the range. */
  if (st_data->force_squarefree) {
    change_lower_strict(t1q, t2q, STATE);
    fmpq_neg(t2q, t2q);
    if (k%2==1) change_upper_strict(t1q, t2q, STATE);
    else change_lower_strict(t1q, t2q, STATE);
  }
  else {
    change_lower(t1q, t2q, STATE);
    fmpq_neg(t2q, t2q);
    if (k%2==1) change_upper(t1q, t2q, STATE);
    else change_lower(t1q, t2q, STATE);
  }
  if (fmpz_cmp(lower, upper) > 0) return(0);

  /* If modulus==0, then return 1 iff [lower, upper] contains 0
     and the Rolle condition is satisfied. */
  if (fmpz_is_zero(modulus)) {
    if ((fmpz_sgn(lower) > 0) || (fmpz_sgn(upper) < 0) ||
	!_fmpz_poly_all_real_roots(POL, NULL)) return(0);
    fmpz_zero(lower);
    fmpz_zero(upper);
    return(1);
  }

  /* Condition: the truncated Hausdorff moment criterion. */
  for (r=0; r<=1; r++) {
    if (k%2 == 1 && !q_is_1) continue; // TODO: remove this restriction
    s = (k%2 == 0 && r == 1) ? 0 : 1;
    fmpq_mat_window_init(hankel_mat, dy_data->hankel_mat, 0, 0, k/2+s, k/2+s);
    for (i=0; i<k/2+s; i++)
      for (j=0; j<k/2+s; j++) {
        t = fmpq_mat_entry(hankel_mat, i, j);
        if (k%2 == 0 && r == 0) 
          fmpq_set(t, fmpq_mat_entry(dy_data->power_sums, i+j, 0));
	else if (k%2 == 0 && r == 1) {
  	  fmpq_mul_ui(t, fmpq_mat_entry(dy_data->power_sums, i+j, 0), 4);
  	  if (!q_is_1) fmpq_mul_fmpz(t, t, q);
	  fmpq_sub(t, t, fmpq_mat_entry(dy_data->power_sums, i+j+2, 0));
	} else {
	  fmpq_mul_ui(t, fmpq_mat_entry(dy_data->power_sums, i+j, 0), 2);
	  if (r == 0) fmpq_add(t, t, fmpq_mat_entry(dy_data->power_sums, i+j+1, 0));
	  else fmpq_sub(t, t, fmpq_mat_entry(dy_data->power_sums, i+j+1, 0));
	}
      } // Final value of t will be used again
    t1 = fmpq_mat_entry(dy_data->hankel_dets[r], k, 0);
    fmpq_mat_det(t1, hankel_mat);
    fmpq_mat_window_clear(hankel_mat);    
    if (k > 1) t2 = fmpq_mat_entry(dy_data->hankel_dets[r], k-2, 0);
    if (k > 1 && fmpq_sgn(t2) > 0) {
      fmpq_div_raw(t1q, t1, t2);
      change_by_sign(t1q, NULL, r, STATE);
    }
    else if ((k > 1 && st_data->force_squarefree) || fmpq_sgn(t1) < 0) return(0);
    else change_by_sign(t, NULL, r, STATE); // t was set in the for loop
    if (fmpz_cmp(lower, upper) > 0) return(0);
  }
  
  /* Condition: tpol has all real roots.
     The range of values where this holds, if nonempty, is an interval. */

  while (1) {
    /* If the range is a singleton, test that value. */
    r = fmpz_cmp(lower, upper);
    if (!r && _fmpz_poly_all_real_roots(POL, lower)) {
      fmpz_set(t2z, lower);
      break;
    }
    if (r >= 0) return(0);

    /* Test the midpoint of the given range. If we find all real roots,
       run a binary search for the left endpoint. */
    fmpz_cmid(t0z, lower, upper);
    if (_fmpz_poly_all_real_roots(POL, t0z)) {
      fmpz_set(t2z, t0z); // Start the right endpoint search from here
      while (fmpz_cmp(lower, t0z)) {
        fmpz_fmid(t1z, lower, t0z);
        if (_fmpz_poly_all_real_roots(POL, t1z)) fmpz_set(t0z, t1z);
        else fmpz_add_ui(lower, t1z, 1); 
      }
      break;
    }

    /* Run a linear search up to the midpoint. */
    do {
      r = _fmpz_poly_all_real_roots(POL, lower);
      if (r) break;
      fmpz_add_ui(lower, lower, 1);
    } while (fmpz_cmp(lower, t0z));
    if (r) {
      fmpz_set(t2z, lower);
      fmpz_sub_ui(upper, t0z, 1); // Use an improved upper bound
      break;
    }
 
    /* Shorten the interval and try again. */
    fmpz_add_ui(lower, t0z, 1);
  }

  /* Use a binary search to find the right endpoint of the interval where the
     real roots condition is satisfied. */
  while (fmpz_cmp(t2z, upper)) {
    fmpz_cmid(t1z, t2z, upper);
    if (_fmpz_poly_all_real_roots(POL, t1z)) fmpz_set(t2z, t1z);
    else fmpz_sub_ui(upper, t1z, 1);
  }

  /* Set the new upper bound. */
  fmpz_mul(upper, upper, modulus);
  fmpz_add(dy_data->upper+n-1, pol+n-1, upper);

  /* Set the new polynomial value, then correct the k-th power sum and related quantities. */
  step_forward(st_data, dy_data, n-1, lower); // see below

  return(1);
}

/* Increment the current moving counter and update stored data to match. */
inline void step_forward(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int n, fmpz_t step) {
  int k = st_data->d-n;
  fmpz *pol = dy_data->pol;
  fmpq *f = st_data->f+n;
  fmpq *t;
  fmpq *t0q = dy_data->w2;
  
  if (step == NULL) {
    fmpq_set(t0q, f);
    fmpz_add(pol+n, pol+n, st_data->modlist+n);
  }
  else {
    fmpq_mul_fmpz(t0q, f, step);
    fmpz_addmul(pol+n, step, st_data->modlist+n);
  }
  t = fmpq_mat_entry(dy_data->power_sums, k, 0);
  fmpq_sub(t, t, t0q);
  t = fmpq_mat_entry(dy_data->hankel_dets[0], k, 0);
  if (k > 1) fmpq_submul(t, fmpq_mat_entry(dy_data->hankel_dets[0], k-2, 0), t0q);
  else fmpq_sub(t, t, t0q);
  t = fmpq_mat_entry(dy_data->hankel_dets[1], k, 0);
  if (k > 1) fmpq_addmul(t, fmpq_mat_entry(dy_data->hankel_dets[1], k-2, 0), t0q);
  else fmpq_add(t, t, t0q);
}

/* Return value sent back in dy_data->flag:
   1: in process
   2: found a solution
   0: tree exhausted
   -1: maximum number of nodes reached
*/

void next_pol(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int max_steps) {
  int d = st_data->d;
  int node_limit = st_data->node_limit;
  fmpz *modlist = st_data->modlist;

  int ascend = dy_data->ascend;
  int n = dy_data->n;
  long node_count = dy_data->node_count;
  fmpz *upper = dy_data->upper;
  fmpz *pol = dy_data->pol;
  fmpz *sympol = dy_data->sympol;
  fmpz *temp = dy_data->w;

  int i, j, flag = 1, count_steps = 0;

  if (dy_data==NULL || !dy_data->flag) return; // No work assigned to this process
  if (n>d) return;

  dy_data->flag = 0; // Prevent work-stealing while this process is running

  while ((flag==1) && (count_steps <= max_steps)) {
    count_steps += 1;
    if (ascend) { // Ascend the tree and step forward as needed.
      n += ascend;
      if (n>d) flag = 0; // This process is complete.
      else {
	ascend = (fmpz_is_zero(modlist+n) || (fmpz_cmp(pol+n, upper+n) >= 0));
	if (!ascend) step_forward(st_data, dy_data, n, NULL);
      }
    } else if (n < 0) { // Return a solution.
      _fmpz_vec_zero(sympol, 2*d+3);
      for (i=0; i<=d; i++) {
	fmpz_set_si(temp, st_data->sign);
	for (j=0; j<=i; j++) {
	  fmpz_addmul(sympol+d+i-2*j, pol+i, temp);
	  if (j<i) {
	    fmpz_mul(temp, temp, st_data->q);
	    fmpz_mul_si(temp, temp, i-j);
	    fmpz_divexact_si(temp, temp, j+1);
	  }
	}
      }
      ascend = 1;
      flag = 2;
    } else { // Compute children of the current node.
      dy_data->n = n;
      ascend = !set_range_from_power_sums(st_data, dy_data);
      n -= 1;
      if (ascend) {
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
