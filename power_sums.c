#include <flint.h>
#include <fmpz_poly.h>
#include <fmpq.h>
#include <fmpq_mat.h>
#include <arith.h>

#include "all_roots_in_interval.h"

typedef struct power_sums_data {
  /* Immutable quantities */
  int d, lead;
  fmpz_t a, b, const1;
  fmpq_mat_t sum_col, sum_prod;
  fmpq_mat_t *sum_mats;
  int *modlist;

  /* Mutable quantities */
  /* Todo: include the polynomial to be tested */

  /* Scratch space */
  fmpz_poly_t pol;  
  fmpz *w;
  slong wlen; /* = 4*d+12 */
  fmpq *w2;
  slong w2len; /* = 7 */
} power_sums_data_t;

void fmpq_floor(fmpz_t res, fmpq_t a) {
  fmpz_fdiv_q(res, fmpq_numref(a), fmpq_denref(a));
};

void fmpq_ceil(fmpz_t res, fmpq_t a) {
  fmpz_cdiv_q(res, fmpq_numref(a), fmpq_denref(a));
};

/* Return the k-th power sum of an integer polynomial as a column vector. */
void fmpq_poly_power_sums(power_sums_data_t *data, int *pol, int k)
{
  int i,j;
  int d = data->d;
  fmpq *t;
  fmpq_t u;

  fmpq_init(u);

  t = fmpq_mat_entry(data->sum_col, 0, 0);
  fmpq_set_si(t, d, 1);
  for (i=1; i<=k; i++) {
    t = fmpq_mat_entry(data->sum_col, i, 0);
    if (i <= d) {
      fmpq_set_si(t, -i*pol[d-i], pol[d]);
    } else {
      fmpq_zero(t);
    }
    for (j=1; j<i && j<=d; j++) {
      fmpq_set_si(u, -pol[d-j], pol[d]);
      fmpq_mul(u, u, fmpq_mat_entry(data->sum_col, 0, i-j));
      fmpq_add(t, t, u);
    }
  }
  for (i=k+1; i<=d; i++) {
    t = fmpq_mat_entry(data->sum_col, i, 0);
    fmpq_zero(t);
  }
  fmpq_clear(u);
}

/* Memory allocation and release.
 */
power_sums_data_t *ranger_init(int d, int lead, int *modlist) {
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
  
  data->modlist = (int *)malloc((d+1)*sizeof(int));
  for (i=0; i<=d; i++) {
    data->modlist[i] = modlist[i];
  }
  
  data->sum_mats = (fmpq_mat_t *)malloc((d+1)*sizeof(fmpq_mat_t));
  data->wlen = 4*d+12;
  data->w = _fmpz_vec_init(data->wlen);
  data->w2len = 7;
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

/* Return values: 1 if bounds[0] <= bounds[1], 0 otherwise. */
int range_from_power_sums(int *bounds, power_sums_data_t *data,
			   int *pol, int k) {
  int i, j, r, modulus, d = data->d;
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
  
  fmpz *tmp1zptr; 
  fmpq *tmp1qptr, *tmp2qptr, *tmp3qptr;

  fmpz *tpol = data->w;
  fmpz *tmp1z = data->w + d + 1;
  
  /* Compute the divided (d-k+1)-st derivative of pol, 
     then determine whether its roots are all in [-2, 2]. */
  fmpz_poly_zero(data->pol);
  for (i=0; i<=k-1; i++) {
    fmpz_set_si(tmp1z, pol[i+d-k+1]);
    for (j=1; j<=d-k+1; j++) {
      fmpz_mul_si(tmp1z, tmp1z, i+j);
      fmpz_divexact_si(tmp1z, tmp1z, j);
    }
    fmpz_poly_set_coeff_fmpz(data->pol, i, tmp1z);
  }

  r = _fmpz_poly_all_roots_in_interval(data->pol->coeffs, k,
				       data->a, data->b, data->w+d+1);

  if (r==0) {
    bounds[0] = 1;
    bounds[1] = 0;
  } else if ((k > d) || (data->modlist[d-k] == 0)) {
    bounds[0] = 0;
    bounds[1] = 0;
  } else {
    modulus = data->modlist[d-k];
    t0z = data->w + d + 2;
    lower = data->w + d + 3;
    upper = data->w + d + 4;
    fmpz *tmp2z = data->w + d + 5;
    fmpz *tmp3z = data->w + d + 6;

    t0q = data->w2;
    f = data->w2 + 1;
    fmpq *tmp1q = data->w2 + 2;
    fmpq *tmp2q = data->w2 + 3;
    fmpq *tmp3q = data->w2 + 4;
    fmpq *tmp4q = data->w2 + 5;
    fmpq *tmp5q = data->w2 + 6;
    
    fmpq_set_si(f, data->lead, modulus*k);
    fmpq_poly_power_sums(data, pol, k);
    fmpq_mat_mul(data->sum_prod, data->sum_mats[k], data->sum_col);

    /* Initial bounds coming from power sums of the asymmetrized polynomial. */
    fmpq_set_si(tmp1q, -2*d, 1);
    fmpq_add(tmp1q, tmp1q, fmpq_mat_entry(data->sum_prod, 0, 0));
    set_lower(tmp1q);
    fmpq_set_si(tmp1q, 2*d, 1);
    fmpq_add(tmp1q, tmp1q, fmpq_mat_entry(data->sum_prod, 0, 0));
    set_upper(tmp1q);

    /* Apply Descartes' rule of signs at data->a = -2, data->b = +2. */
    /* Currently data->pol is the divided (d-k+1)-st derivative of pol. */
    fmpz_poly_evaluate_fmpz(tmp3z, data->pol, data->a);
    fmpq_set_fmpz_frac(tmp3q, tmp3z, data->const1);
    fmpz_poly_evaluate_fmpz(tmp3z, data->pol, data->b);
    fmpq_set_fmpz_frac(tmp4q, tmp3z, data->const1);

    /* Set data->pol to the divided (d-k)-th derivative of pol. */
    fmpz_poly_shift_left(data->pol, data->pol, 1);
    fmpz_poly_set_coeff_ui(data->pol, 0, pol[d-k]);
    for (i=1; i<=k; i++) {
      tmp1zptr = fmpz_poly_get_coeff_ptr(data->pol, i);
      fmpz_mul_si(tmp1zptr, tmp1zptr, d-k+1);
      fmpz_divexact_si(tmp1zptr, tmp1zptr, i);
    }

    fmpz_poly_evaluate_fmpz(tmp3z, data->pol, data->a);
    fmpq_set_fmpz_frac(tmp1q, tmp3z, data->const1);
    fmpz_poly_evaluate_fmpz(tmp3z, data->pol, data->b);
    fmpq_set_fmpz_frac(tmp2q, tmp3z, data->const1);

    fmpq_set_si(tmp5q, -k, data->lead);
    fmpq_mul(tmp1q, tmp1q, tmp5q);
    r = fmpq_sgn(tmp3q);
    if (r >= 0) {
      change_upper(tmp1q);
    }
    if (r <= 0) {
      change_lower(tmp1q);
    }
    fmpq_mul(tmp2q, tmp2q, tmp5q);
    r = fmpq_sgn(tmp4q);
    if (r >= 0) {
      change_lower(tmp2q);
    }
    if (r <= 0) {
      change_upper(tmp2q);
    }
    
    /* Additional bounds based on power sums. */
    if ((fmpz_cmp(lower, upper) <= 0) && k >= 2) {
      fmpq_add(tmp1q, fmpq_mat_entry(data->sum_prod, 1, 0),
	       fmpq_mat_entry(data->sum_prod, 2, 0));
      fmpq_set_si(tmp2q, 4*d, 1);
      fmpq_sub(tmp3q, tmp1q, tmp2q);
      if (k==2) {
	fmpq_set_si(tmp4q, 1, 2);
	fmpq_mul(tmp3q, tmp3q, tmp4q);
      }
      change_lower(tmp3q);
      fmpq_add(tmp3q, tmp1q, tmp2q);
      if (k==2) {
	fmpq_mul(tmp3q, tmp3q, tmp4q);
      }
      change_upper(tmp3q);

      tmp1qptr = fmpq_mat_entry(data->sum_prod, 3, 0);
      tmp2qptr = fmpq_mat_entry(data->sum_prod, 4, 0);
      tmp3qptr = fmpq_mat_entry(data->sum_prod, 5, 0);
      if (fmpq_sgn(tmp3qptr) > 0) {
	fmpq_mul(tmp4q, tmp2qptr, tmp2qptr);
	fmpq_div(tmp4q, tmp4q, tmp3qptr);
	fmpq_sub(tmp4q, tmp1qptr, tmp4q);
	change_upper(tmp4q);
      }
      fmpq_set_si(tmp4q, -4, 1);
      fmpq_mul(tmp4q, tmp4q, tmp2qptr);
      fmpq_add(tmp4q, tmp4q, tmp1qptr);
      change_lower(tmp4q);
	
      tmp1qptr = fmpq_mat_entry(data->sum_prod, 6, 0);
      tmp2qptr = fmpq_mat_entry(data->sum_prod, 7, 0);
      tmp3qptr = fmpq_mat_entry(data->sum_prod, 8, 0);
      fmpq_mul(tmp4q, tmp2qptr, tmp2qptr);
      if ((k%2 == 0) && (fmpq_sgn(tmp3qptr) > 0)) {
	fmpq_div(tmp4q, tmp4q, tmp3qptr);
	fmpq_sub(tmp4q, tmp1qptr, tmp4q);
	change_upper(tmp4q);
      } else if ((k%2 == 1) && (fmpq_sgn(tmp3qptr) < 0)) {
	fmpq_div(tmp4q, tmp4q, tmp3qptr);
	fmpq_sub(tmp4q, tmp1qptr, tmp4q);
	change_lower(tmp4q);
      }
      fmpq_set_si(tmp4q, 4, 1);
      fmpq_mul(tmp4q, tmp4q, tmp2qptr);
      fmpq_add(tmp4q, tmp4q, tmp1qptr);
      if (k%2 == 0) {
	change_lower(tmp4q);
      } else {
	change_upper(tmp4q);
      }

      if (k%2 == 0) {
	fmpq_set_si(tmp4q, -4, 1);
	fmpq_mul(tmp4q, tmp4q, fmpq_mat_entry(data->sum_col, k-2, 0));
	fmpq_add(tmp4q, tmp4q, fmpq_mat_entry(data->sum_col, k, 0));
	change_lower(tmp4q);	
      }
    }
    
    bounds[0] = fmpz_get_si(lower)*modulus;
    bounds[1] = fmpz_get_si(upper)*modulus;
  }
  return(bounds[0] <= bounds[1]);

}
