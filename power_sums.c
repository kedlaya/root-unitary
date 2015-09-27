#include <flint.h>
#include <fmpz_poly.h>
#include <fmpq.h>
#include <fmpq_mat.h>
#include <arith.h>

#include "all_roots_in_interval.h"

typedef struct power_sums_data {
  int d, lead;
  fmpq_mat_t *sum_mats;
  fmpq_mat_t sum_col, sum_prod;
  fmpz *w;
  slong wlen;
} power_sums_data_t;

void fmpq_floor(fmpz_t res, fmpq_t a) {
  fmpz_fdiv_q(res, fmpq_numref(a), fmpq_denref(a));
};

void fmpq_ceil(fmpz_t res, fmpq_t a) {
  fmpz_cdiv_q(res, fmpq_numref(a), fmpq_denref(a));
};

/* For some reason, I can't import this from FLINT. */
void chebyshev_t(fmpz_poly_t pol, int n) {
  int k;
  fmpz_t m;
  fmpz_poly_t u1, u2;

  fmpz_init(m);
  fmpz_poly_init(u1);
  fmpz_poly_init(u2);
  
  fmpz_poly_zero(pol);

  fmpz_poly_zero(u1);
  fmpz_poly_set_coeff_si(u1, 0, -1);
  fmpz_poly_set_coeff_si(u1, 2, 1);
  
  for (k = 0; k <= n/2; k++) {
    fmpz_poly_pow(u2, u1, k);
    fmpz_poly_shift_left(u2, u2, n-2*k);
    fmpz_bin_uiui(m, n, 2*k);
    fmpz_poly_scalar_mul_fmpz(u2, u2, m);
    fmpz_poly_add(pol, pol, u2);
  }

  fmpz_clear(m);
  fmpz_poly_clear(u1);
  fmpz_poly_clear(u2);
}

/* Return the k-th power sum of an integer polynomial as a column vector. */
void fmpq_poly_power_sums(power_sums_data_t data, int *pol, int k)
{
  int i,j;
  fmpq *t;
  fmpq_t u;

  fmpq_init(u);

  t = fmpq_mat_entry(data.sum_col, 0, 0);
  fmpq_set_si(t, data.d, 1);
  for (i=1; i<=k; i++) {
    t = fmpq_mat_entry(data.sum_col, i, 0);
    if (i <= data.d) {
      fmpq_set_si(t, -i*pol[data.d-i], pol[data.d]);
    } else {
      fmpq_zero(t);
    }
    for (j=1; j<i && j<=data.d; j++) {
      fmpq_set_si(u, -pol[data.d-j], pol[data.d]);
      fmpq_mul(u, u, fmpq_mat_entry(data.sum_col, 0, i-j));
      fmpq_add(t, t, u);
    }
  }
  for (i=k+1; i<=data.d; i++) {
    t = fmpq_mat_entry(data.sum_col, i, 0);
    fmpq_zero(t);
  }
  fmpq_clear(u);
}

/* Memory allocation and release.
 */
power_sums_data_t ranger_init(int d, int lead, int *modlist) {
  int i, j;
  power_sums_data_t data;
  fmpz_poly_t pol;
  fmpz_t l, m;
  fmpq *k1;
  
  fmpz_poly_init(pol);
  fmpz_init(l);
  fmpz_init(m);

  fmpz_one(l);  
  
  data.d = d;
  data.lead = lead;
  data.sum_mats = (fmpq_mat_t *)malloc((d+1)*sizeof(fmpq_mat_t));
  fmpq_mat_init(data.sum_col, d+1, 1);
  fmpq_mat_init(data.sum_prod, 9, 1);
  for (i=0; i<=d; i++) {

    fmpq_mat_init(data.sum_mats[i], 9, d+1);

    /* Compute Chebyshev polynomial of the first kind. */
    chebyshev_t(pol, i);
    for (j=0; j<=d; j++) {
      /* Row 0: coeffs of 2*(i-th Chebyshev polynomial)(x/2). */
      k1 = fmpq_mat_entry(data.sum_mats[i], 0, j);
      if (j > i) {
	fmpq_zero(k1);
      }
      else {
	fmpq_set_fmpz_frac(k1, fmpz_poly_get_coeff_ptr(pol, j), l);
	fmpz_mul_2exp(m, l, j);
	fmpq_div_fmpz(k1, k1, m);
	fmpz_set_ui(m, 2);
	fmpq_mul_fmpz(k1, k1, m);
      }

      /* Row 1: coeffs of row 0 from matrix i-2, multiplied by -2. */
      k1 = fmpq_mat_entry(data.sum_mats[i], 1, j);
      if (i < 2) {
	fmpq_zero(k1);
      } else {
	fmpq_set(k1, fmpq_mat_entry(data.sum_mats[i-2], 0, j));
	fmpz_set_si(m, -2);
	fmpq_mul_fmpz(k1, k1, m);
      }

      /* Row 2: coeffs of row 0 from matrix i-2, shifted by 2. */
      k1 = fmpq_mat_entry(data.sum_mats[i], 2, j);
      if (i<2 || j<2) {
	fmpq_zero(k1);
      } else {
	fmpq_set(k1, fmpq_mat_entry(data.sum_mats[i-2], 0, j-2));
      }

      /* Row 3: coeffs of (2+x)^i. */
      k1 = fmpq_mat_entry(data.sum_mats[i], 3, j);
      if (j>i) {
	fmpq_zero(k1);
      } else {
	fmpz_bin_uiui(m, i, j);
	fmpq_set_fmpz_frac(k1, m, l);
	fmpq_mul_2exp(k1, k1, i-j);
      }
      
      /* Row 4: coeffs of (2+x)^(i-1). */
      k1 = fmpq_mat_entry(data.sum_mats[i], 4, j);
      if (i<1) {
	fmpq_zero(k1);
      } else {
	fmpq_set(k1, fmpq_mat_entry(data.sum_mats[i-1], 3, j));	
      }

      /* Row 5: coeffs of (2+x)^(i-2). */
      k1 = fmpq_mat_entry(data.sum_mats[i], 5, j);
      if (i<2) {
	fmpq_zero(k1);
      } else {
	fmpq_set(k1, fmpq_mat_entry(data.sum_mats[i-2], 3, j));	
      }

      /* Row 6: coeffs of (-2+x)^i. */
      k1 = fmpq_mat_entry(data.sum_mats[i], 6, j);
      fmpq_set(k1, fmpq_mat_entry(data.sum_mats[i], 3, j));
      if ((i-j)%2==1) {
	fmpq_neg(k1, k1);
      }

      /* Row 7: coeffs of (-2+x)^(i-1). */
      k1 = fmpq_mat_entry(data.sum_mats[i], 7, j);
      if (i<1) {
	fmpq_zero(k1);
      } else {
	fmpq_set(k1, fmpq_mat_entry(data.sum_mats[i-1], 6, j));	
      }

      /* Row 8: coeffs of (-2+x)^(i-2). */
      k1 = fmpq_mat_entry(data.sum_mats[i], 8, j);
      if (i<2) {
	fmpq_zero(k1);
      } else {
	fmpq_set(k1, fmpq_mat_entry(data.sum_mats[i-2], 6, j));	
      }
    }
  }
  
  fmpz_poly_clear(pol);
  fmpz_clear(l);
  fmpz_clear(m);
  return(data);
}

void ranger_clear(power_sums_data_t data) {
  int i;
  for (i=0; i<=data.d; i++) {
    fmpq_mat_clear(data.sum_mats[i]);
  }
  free(data.sum_mats);
  fmpq_mat_clear(data.sum_col);
}

/* Return values: 1 if bounds[0] <= bounds[1], 0 otherwise. */
int range_from_power_sums(int *bounds, power_sums_data_t data,
			   int *pol, int modulus, int k) {
  int i, j, r, d = data.d;
  fmpz_t lower, upper, tmp0z;
  fmpq_t tmp0q, f; 
    
  void set_lower(fmpq_t val) {
    fmpq_mul(tmp0q, val, f);
    fmpq_ceil(lower, tmp0q);
  }
  
  void set_upper(fmpq_t val) {
    fmpq_mul(tmp0q, val, f);
    fmpq_floor(upper, tmp0q);
  }

  void change_lower(fmpq_t val) {
    fmpq_mul(tmp0q, val, f);
    fmpq_ceil(tmp0z, tmp0q);
    if (fmpz_cmp(tmp0z, lower) > 0) {
      fmpz_set(lower, tmp0z);
    }
  }
  
  void change_upper(fmpq_t val) {
    fmpq_mul(tmp0q, val, f);
    fmpq_floor(tmp0z, tmp0q);
    if (fmpz_cmp(tmp0z, upper) < 0) {
      fmpz_set(upper, tmp0z);
    }
  }
  
  fmpz_t tmp1z, tmp2z, tmp3z, tmp4z;
  fmpq_t tmp1q, tmp2q, tmp3q, tmp4q, tmp5q;

  fmpz_poly_t tmpzpol;
  fmpq *tmp1qptr, *tmp2qptr, *tmp3qptr;
  fmpz *w;

  fmpz_init(tmp1z);
  fmpz_init(tmp2z);
  fmpz_poly_init(tmpzpol);

  /* Compute the divided (d-k+1)-st derivative of pol, 
     then determine whether its roots are all in [-2, 2]. */
  for (i=0; i<=k-1; i++) {
    fmpz_set_si(tmp1z, pol[i+d-k+1]);
    for (j=1; j<=d-k+1; j++) {
      fmpz_mul_si(tmp1z, tmp1z, i+j);
      fmpz_divexact_si(tmp1z, tmp1z, j);
    }
    fmpz_poly_set_coeff_fmpz(tmpzpol, i, tmp1z);
  }
  fmpz_set_si(tmp1z, -2);
  fmpz_set_si(tmp2z, 2);
  const slong n = tmpzpol->length;
  w = _fmpz_vec_init(3 * n + 9);
  r = _fmpz_poly_all_roots_in_interval(tmpzpol->coeffs, n, tmp1z, tmp2z, w);
  _fmpz_vec_clear(w, 3 * n + 9);

  if (r==0) {
    bounds[0] = 1;
    bounds[1] = 0;
  } else if ((k > d) || (modulus == 0)) {
    bounds[0] = 0;
    bounds[1] = 0;
  } else {
    fmpz_init(lower);
    fmpz_init(upper);
    fmpz_init(tmp0z);
    fmpz_init(tmp3z);
    fmpz_init(tmp4z);

    fmpq_init(tmp0q);
    fmpq_init(f);
    fmpq_init(tmp1q);
    fmpq_init(tmp2q);
    fmpq_init(tmp3q);
    fmpq_init(tmp4q);
    fmpq_init(tmp5q);
    
    fmpq_set_si(f, data.lead, modulus*k);
    fmpq_poly_power_sums(data, pol, k);
    fmpq_mat_mul(data.sum_prod, data.sum_mats[k], data.sum_col);

    /* Initial bounds coming from power sums of the asymmetrized polynomial. */
    fmpq_set_si(tmp1q, -2*d, 1);
    fmpq_add(tmp1q, tmp1q, fmpq_mat_entry(data.sum_prod, 0, 0));
    set_lower(tmp1q);
    fmpq_set_si(tmp1q, 2*d, 1);
    fmpq_add(tmp1q, tmp1q, fmpq_mat_entry(data.sum_prod, 0, 0));
    set_upper(tmp1q);

    /* Apply Descartes' rule of signs at tmp1z = -2, tmp2z = +2. */
    /* Currently tmpzpol is the divided (d-k+1)-st derivative of pol. */
    fmpz_one(tmp4z);
    fmpz_poly_evaluate_fmpz(tmp3z, tmpzpol, tmp1z);
    fmpq_set_fmpz_frac(tmp3q, tmp3z, tmp4z);
    fmpz_poly_evaluate_fmpz(tmp3z, tmpzpol, tmp2z);
    fmpq_set_fmpz_frac(tmp4q, tmp3z, tmp4z);

    /* Set tmpzpol to the divided (d-k)-th derivative of pol. */
    fmpz_poly_shift_left(tmpzpol, tmpzpol, 1);
    fmpz_poly_set_coeff_ui(tmpzpol, 0, pol[d-k]);
    for (i=1; i<=k; i++) { // Todo: manipulate coefficient by reference.
      fmpz_poly_get_coeff_fmpz(tmp3z, tmpzpol, i);
      fmpz_mul_si(tmp3z, tmp3z, d-k+1);
      fmpz_divexact_si(tmp3z, tmp3z, i);
      fmpz_poly_set_coeff_fmpz(tmpzpol, i, tmp3z);
    }
    fmpz_poly_evaluate_fmpz(tmp3z, tmpzpol, tmp1z);
    fmpq_set_fmpz_frac(tmp1q, tmp3z, tmp4z);
    fmpz_poly_evaluate_fmpz(tmp3z, tmpzpol, tmp2z);
    fmpq_set_fmpz_frac(tmp2q, tmp3z, tmp4z);

    fmpq_set_si(tmp5q, -k, data.lead);
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
      fmpq_add(tmp1q, fmpq_mat_entry(data.sum_prod, 1, 0),
	       fmpq_mat_entry(data.sum_prod, 2, 0));
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

      tmp1qptr = fmpq_mat_entry(data.sum_prod, 3, 0);
      tmp2qptr = fmpq_mat_entry(data.sum_prod, 4, 0);
      tmp3qptr = fmpq_mat_entry(data.sum_prod, 5, 0);
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
	
      tmp1qptr = fmpq_mat_entry(data.sum_prod, 6, 0);
      tmp2qptr = fmpq_mat_entry(data.sum_prod, 7, 0);
      tmp3qptr = fmpq_mat_entry(data.sum_prod, 8, 0);
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
	fmpq_mul(tmp4q, tmp4q, fmpq_mat_entry(data.sum_col, k-2, 0));
	fmpq_add(tmp4q, tmp4q, fmpq_mat_entry(data.sum_col, k, 0));
	change_lower(tmp4q);	
      }
    }
    
    bounds[0] = fmpz_get_si(lower)*modulus;
    bounds[1] = fmpz_get_si(upper)*modulus;
    
    fmpz_clear(lower);
    fmpz_clear(upper);
    fmpz_clear(tmp0z);
    fmpz_clear(tmp3z);
    fmpz_clear(tmp4z);

    fmpq_clear(tmp0q);
    fmpq_clear(f);
    fmpq_clear(tmp1q);
    fmpq_clear(tmp2q);
    fmpq_clear(tmp3q);
    fmpq_clear(tmp4q);
    fmpq_clear(tmp5q);
    
  }
  fmpz_clear(tmp1z);
  fmpz_clear(tmp2z);
  fmpz_poly_clear(tmpzpol);
  return(bounds[0] <= bounds[1]);

}
