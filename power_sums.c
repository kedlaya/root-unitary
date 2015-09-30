#include <flint.h>
#include <fmpz_poly.h>
#include <fmpq.h>
#include <fmpq_mat.h>
#include <arith.h>

#include "all_roots_in_interval.h"

typedef struct power_sums_data {
  /* Immutable quantities */
  int d, lead;
  int *modlist;
  fmpz_t a, b, const1;
  fmpz_mat_t binom_mat;
  fmpq_mat_t *sum_mats;

  /* Mutable quantities */
  /* Todo: include the polynomial to be tested */
  fmpq_mat_t sum_col, sum_prod;
  fmpz_poly_t pol;  

  /* Scratch space */
  fmpz *w;
  slong wlen; /* = 4*d+12 */
  fmpq *w2;
  slong w2len; /* = 5 */
} power_sums_data_t;

void fmpq_floor(fmpz_t res, fmpq_t a) {
  fmpz_fdiv_q(res, fmpq_numref(a), fmpq_denref(a));
};

void fmpq_ceil(fmpz_t res, fmpq_t a) {
  fmpz_cdiv_q(res, fmpq_numref(a), fmpq_denref(a));
};

/* Memory allocation and release.
 */
power_sums_data_t *ranger_init(int d, int lead, int *Q0, int *modlist) {
  int i, j;
  power_sums_data_t *data;
  fmpz_poly_t pol;
  fmpz_t m;
  fmpq *k1;

  fmpz_poly_init(pol);
  fmpz_init(m);

  data = (power_sums_data_t *)malloc(sizeof(power_sums_data_t));
  data->d = d;
  data->lead = lead;

  fmpz_init(data->a);
  fmpz_init(data->b);
  fmpz_init(data->const1);
  fmpz_set_si(data->a, -2);
  fmpz_set_si(data->b, 2);
  fmpz_one(data->const1);
  fmpz_poly_init(data->pol);
  for (i=0; i<=d; i++) {
    fmpz_poly_set_coeff_si(data->pol, i, Q0[i]);
  }
  
  data->modlist = (int *)malloc((d+1)*sizeof(int));
  for (i=0; i<=d; i++) {
    data->modlist[i] = modlist[i];
  }
  fmpz_mat_init(data->binom_mat, d+1, d+1);
  for (i=0; i<=d; i++) {
    for (j=0; j<=d; j++) {
      fmpz_bin_uiui(fmpq_mat_entry(data->binom_mat, i, j), i, j);
    }
  }
  
  data->sum_mats = (fmpq_mat_t *)malloc((d+1)*sizeof(fmpq_mat_t));
  data->wlen = 4*d+12;
  data->w = _fmpz_vec_init(data->wlen);
  data->w2len = 5;
  data->w2 = _fmpq_vec_init(data->w2len);
  fmpq_mat_init(data->sum_col, d+1, 1);
  fmpq_mat_init(data->sum_prod, 9, 1);
  for (i=0; i<=d; i++) {

    fmpq_mat_init(data->sum_mats[i], 9, d+1);

    arith_chebyshev_t_polynomial(pol, i);
    for (j=0; j<=d; j++) {
      
      /* Row 0: coeffs of 2*(i-th Chebyshev polynomial)(x/2). */
      k1 = fmpq_mat_entry(data->sum_mats[i], 0, j);
      if (j > i) {
	fmpq_zero(k1);
      }
      else {
	fmpq_set_fmpz_frac(k1, fmpz_poly_get_coeff_ptr(pol, j), data->const1);
	fmpz_mul_2exp(m, data->const1, j);
	fmpq_div_fmpz(k1, k1, m);
	fmpz_set_ui(m, 2);
	fmpq_mul_fmpz(k1, k1, m);
      }

      /* Row 1: coeffs of row 0 from matrix i-2, multiplied by -2. */
      k1 = fmpq_mat_entry(data->sum_mats[i], 1, j);
      if (i < 2) {
	fmpq_zero(k1);
      } else {
	fmpq_set(k1, fmpq_mat_entry(data->sum_mats[i-2], 0, j));
	fmpz_set_si(m, -2);
	fmpq_mul_fmpz(k1, k1, m);
      }

      /* Row 2: coeffs of row 0 from matrix i-2, shifted by 2. */
      k1 = fmpq_mat_entry(data->sum_mats[i], 2, j);
      if (i<2 || j<2) {
	fmpq_zero(k1);
      } else {
	fmpq_set(k1, fmpq_mat_entry(data->sum_mats[i-2], 0, j-2));
      }

      /* Row 3: coeffs of (2+x)^i. */
      k1 = fmpq_mat_entry(data->sum_mats[i], 3, j);
      if (j>i) {
	fmpq_zero(k1);
      } else {
	fmpz_bin_uiui(m, i, j);
	fmpq_set_fmpz_frac(k1, m, data->const1);
	fmpq_mul_2exp(k1, k1, i-j);
      }
      
      /* Row 4: coeffs of (2+x)^(i-1). */
      k1 = fmpq_mat_entry(data->sum_mats[i], 4, j);
      if (i<1) {
	fmpq_zero(k1);
      } else {
	fmpq_set(k1, fmpq_mat_entry(data->sum_mats[i-1], 3, j));	
      }

      /* Row 5: coeffs of (2+x)^(i-2). */
      k1 = fmpq_mat_entry(data->sum_mats[i], 5, j);
      if (i<2) {
	fmpq_zero(k1);
      } else {
	fmpq_set(k1, fmpq_mat_entry(data->sum_mats[i-2], 3, j));	
      }

      /* Row 6: coeffs of (-2+x)^i. */
      k1 = fmpq_mat_entry(data->sum_mats[i], 6, j);
      fmpq_set(k1, fmpq_mat_entry(data->sum_mats[i], 3, j));
      if ((i-j)%2==1) {
	fmpq_neg(k1, k1);
      }

      /* Row 7: coeffs of (-2+x)^(i-1). */
      k1 = fmpq_mat_entry(data->sum_mats[i], 7, j);
      if (i<1) {
	fmpq_zero(k1);
      } else {
	fmpq_set(k1, fmpq_mat_entry(data->sum_mats[i-1], 6, j));	
      }

      /* Row 8: coeffs of (-2+x)^(i-2). */
      k1 = fmpq_mat_entry(data->sum_mats[i], 8, j);
      if (i<2) {
	fmpq_zero(k1);
      } else {
	fmpq_set(k1, fmpq_mat_entry(data->sum_mats[i-2], 6, j));	
      }
    }
  }
  
  fmpz_poly_clear(pol);
  fmpz_clear(m);
  return(data);
}

void ranger_clear(power_sums_data_t *data) {
  int i;
  fmpz_clear(data->a);
  fmpz_clear(data->b);
  fmpz_init(data->const1);
  fmpz_poly_clear(data->pol);
  fmpz_mat_clear(data->binom_mat);
  for (i=0; i<=data->d; i++) {
    fmpq_mat_clear(data->sum_mats[i]);
  }
  free(data->sum_mats);
  fmpq_mat_clear(data->sum_col);
  free(data->modlist);
  _fmpz_vec_clear(data->w, data->wlen);
  _fmpq_vec_clear(data->w2, data->w2len);
  free(data);
}

/* Compute the 0th to k-th power sums and place them in data->sum_col. */
void fmpq_poly_power_sums(power_sums_data_t *data, int *pol, int k)
{
  int i,j;
  int d = data->d;
  fmpq *t, *u;

  u = data->w + d + 1;
  fmpq_set_si(fmpq_mat_entry(data->sum_col, 0, 0), d, 1);
  for (i=1; i<=k; i++) {
    t = fmpq_mat_entry(data->sum_col, i, 0);
    if (i <= d) fmpq_set_si(t, -i*pol[d-i], pol[d]);
    else fmpq_zero(t);
    for (j=1; j<i && j<=d; j++) {
      fmpq_set_si(u, -pol[d-j], pol[d]);
      fmpq_addmul(t, u, fmpq_mat_entry(data->sum_col, 0, i-j));
    }
  }
  for (i=k+1; i<=d; i++) 
    fmpq_zero(fmpq_mat_entry(data->sum_col, i, 0));
}

/* Return values: 1 if bounds[0] <= bounds[1], 0 otherwise. */
int range_from_power_sums(int *bounds, power_sums_data_t *data,
			   int *pol, int k) {
  int i, j, r, r1, r2, modulus, d = data->d;
  fmpz *lower, *upper, *t0z;
  fmpq *t0q, *f;
    
  void set_lower(const fmpq_t val) {
    fmpq_mul(t0q, val, f);
    fmpq_ceil(lower, t0q);
  }
  
  void set_upper(const fmpq_t val) {
    fmpq_mul(t0q, val, f);
    fmpq_floor(upper, t0q);
  }

  void change_lower(const fmpq_t val) {
    fmpq_mul(t0q, val, f);
    fmpq_ceil(t0z, t0q);
    if (fmpz_cmp(t0z, lower) > 0) {
      fmpz_set(lower, t0z);
    }
  }
  
  void change_upper(const fmpq_t val) {
    fmpq_mul(t0q, val, f);
    fmpq_floor(t0z, t0q);
    if (fmpz_cmp(t0z, upper) < 0) {
      fmpz_set(upper, t0z);
    }
  }
  
  fmpz *tpol = data->w;
  t0z = data->w + d + 2;
  
  /* Compute the divided (d-k+1)-st derivative of pol, 
     then determine whether its roots are all in [-2, 2]. */
  for (i=0; i<=k-1; i++)
    fmpz_mul_si(tpol+i, fmpz_mat_entry(data->binom_mat, i+d-k+1, d-k+1), 
		pol[i+d-k+1]);

  r = _fmpz_poly_all_roots_in_interval(tpol, k, data->a, data->b, data->w+d+1);

  /* If r=0, abort. 
     If r=1 and k>d, return trivially; there are no further coefficients to find.
     If r=1 and k<=d, continue to compute bounds.
  */

  if (r==0) {
    bounds[0] = 1;
    bounds[1] = 0;
  } else if ((k > d) || (data->modlist[d-k] == 0)) {
    bounds[0] = 0;
    bounds[1] = 0;
  } else {
    /* Allocate temporary variables from persistent scratch space. */
    lower = data->w + d + 3;
    upper = data->w + d + 4;

    t0q = data->w2;
    f = data->w2 + 1;
    fmpq *t1q = data->w2 + 2;
    fmpq *t2q = data->w2 + 3;
    fmpq *t3q = data->w2 + 4;
    
    /* Initialize bounds using power sums of the asymmetrized polynomial. */
    fmpq_poly_power_sums(data, pol, k);
    modulus = data->modlist[d-k];
    fmpq_set_si(f, data->lead, modulus*k);
    fmpq_mat_mul(data->sum_prod, data->sum_mats[k], data->sum_col);

    fmpq_set_si(t1q, 2*d, 1);
    fmpq_sub(t0q, fmpq_mat_entry(data->sum_prod, 0, 0), t1q);
    set_lower(t0q);
    fmpq_add(t0q, fmpq_mat_entry(data->sum_prod, 0, 0), t1q);
    set_upper(t0q);

    /* Apply Descartes' rule of signs at data->a = -2, data->b = +2. */
    /* Currently tpol is the divided (d-k+1)-st derivative of pol;
       evaluate at the endpoints. */
    _fmpz_poly_evaluate_fmpz(t0z, tpol, k, data->a);
    r1 = fmpz_sgn(t0z);
    _fmpz_poly_evaluate_fmpz(t0z, tpol, k, data->b);
    r2 = fmpz_sgn(t0z);

    /* Now set tpol to the divided (d-k)-th derivative of pol.
       then evaluate again. */
    for (i=k; i>=1; i--) {
      fmpz_set(tpol+i, tpol+i-1);
      fmpz_mul_si(tpol+i, tpol+i, d-k+1);
      fmpz_divexact_si(tpol+i, tpol+i, i);      
    }
    fmpq_set_si(t2q, -k, data->lead);

    fmpz_set_si(tpol, pol[d-k]);    
    _fmpz_poly_evaluate_fmpz(t0z, tpol, k+1, data->a);
    fmpq_mul_fmpz(t1q, t2q, t0z);
    if (r1 >= 0) change_upper(t1q);
    if (r1 <= 0) change_lower(t1q);

    _fmpz_poly_evaluate_fmpz(t0z, tpol, k+1, data->b);
    fmpq_mul_fmpz(t1q, t2q, t0z);
    if (r2 >= 0) change_lower(t1q);
    if (r2 <= 0) change_upper(t1q);
    
    /* Additional bounds based on power sums. */
    if ((fmpz_cmp(lower, upper) <= 0) && k >= 2) {
      fmpq_add(t1q, fmpq_mat_entry(data->sum_prod, 1, 0),
	       fmpq_mat_entry(data->sum_prod, 2, 0));
      fmpq_set_si(t2q, 4*d, 1);
      fmpq_sub(t0q, t1q, t2q);
      // fmpq_sub_si(t3q, t1q, 4*d);
      if (k==2) fmpq_div_fmpz(t0q, t0q, data->b); // b=2
      change_lower(t0q);
      fmpq_add(t0q, t1q, t2q);
      // fmpq_add_si(t3q, t1q, 4*d);
      if (k==2) fmpq_div_fmpz(t0q, t0q, data->b); // b=2
      change_upper(t0q);

      /* t1q, t2q, t3q are no longer needed, so can be reassigned. */
      t1q = fmpq_mat_entry(data->sum_prod, 3, 0);
      t2q = fmpq_mat_entry(data->sum_prod, 4, 0);
      t3q = fmpq_mat_entry(data->sum_prod, 5, 0);
      if (fmpq_sgn(t3q) > 0) { // t0q <- t1q - t2q^2/t3q
	fmpq_mul(t0q, t2q, t2q);
	fmpq_div(t0q, t0q, t3q);
	fmpq_sub(t0q, t1q, t0q);
	change_upper(t0q);
      }
      fmpq_set_si(t3q, -4, 1);
      fmpq_mul(t0q, t3q, t2q);
      // fmpq_mul_si(t0q, t2q, -4);
      fmpq_add(t0q, t0q, t1q);
      change_lower(t0q);
	
      t1q = fmpq_mat_entry(data->sum_prod, 6, 0);
      t2q = fmpq_mat_entry(data->sum_prod, 7, 0);
      t3q = fmpq_mat_entry(data->sum_prod, 8, 0);
      fmpq_mul(t0q, t2q, t2q);
      if ((k%2 == 0) && (fmpq_sgn(t3q) > 0)) {
	fmpq_div(t0q, t0q, t3q);
	fmpq_sub(t0q, t1q, t0q);
	change_upper(t0q);
      } else if ((k%2 == 1) && (fmpq_sgn(t3q) < 0)) {
	fmpq_div(t0q, t0q, t3q);
	fmpq_sub(t0q, t1q, t0q);
	change_lower(t0q);
      }
      fmpq_set_si(t0q, 4, 1);
      fmpq_mul(t0q, t0q, t2q);
      fmpq_add(t0q, t0q, t1q);
      if (k%2 == 0) {
	change_lower(t0q);
      } else {
	change_upper(t0q);
      }

      if (k%2 == 0) {
	fmpq_set_si(t0q, -4, 1);
	fmpq_mul(t0q, t0q, fmpq_mat_entry(data->sum_col, k-2, 0));
	fmpq_add(t0q, t0q, fmpq_mat_entry(data->sum_col, k, 0));
	change_lower(t0q);	
      }
    }
    
    bounds[0] = fmpz_get_si(lower)*modulus;
    bounds[1] = fmpz_get_si(upper)*modulus;
  }
  return(bounds[0] <= bounds[1]);

}
