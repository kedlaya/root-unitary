#include <flint.h>
#include <fmpz_poly.h>
#include <fmpq.h>
#include <fmpq_mat.h>
#include <arith.h>

#include "all_roots_in_interval.h"

typedef struct ps_static_data {
  int d, lead, verbosity, node_count;
  fmpz_t a, b;
  fmpz_mat_t binom_mat;
  fmpz *modlist;
  fmpq_mat_t *sum_mats;
  fmpq_t *f;
} ps_static_data_t;

typedef struct ps_dynamic_data {
  int d, n, count, ascend;
  fmpq_mat_t sum_col, sum_prod;
  fmpz *pol, *upper;

  /* Scratch space */
  fmpz *w;
  slong wlen; /* = 4*d+12 */
  fmpq *w2;
  slong w2len; /* = 5 */
} ps_dynamic_data_t;

void fmpq_floor(fmpz_t res, fmpq_t a) {
  fmpz_fdiv_q(res, fmpq_numref(a), fmpq_denref(a));
};

void fmpq_ceil(fmpz_t res, fmpq_t a) {
  fmpz_cdiv_q(res, fmpq_numref(a), fmpq_denref(a));
};

/* Memory allocation and release.
 */
ps_static_data_t *ps_static_init(int d, int lead, int *modlist, 
				 int verbosity, int node_count) {
  int i, j;
  ps_static_data_t *st_data;
  fmpz_poly_t pol;
  fmpz_t m, const1;
  fmpq *k1;

  fmpz_poly_init(pol);
  fmpz_init(m);
  fmpz_init_set_ui(const1, 1);

  st_data = (ps_static_data_t *)malloc(sizeof(ps_static_data_t));

  st_data->d = d;
  st_data->lead = lead;
  st_data->verbosity = verbosity;
  st_data->node_count = node_count;

  fmpz_init(st_data->a);
  fmpz_init(st_data->b);
  fmpz_set_si(st_data->a, -2);
  fmpz_set_si(st_data->b, 2);

  st_data->modlist =_fmpz_vec_init(d+1);
  st_data->f = _fmpq_vec_init(d+1);
  for (i=0; i<=d; i++) {
    fmpz_set_si(st_data->modlist+i, modlist[i]);
    fmpq_set_si(st_data->f+i, d-i, lead);
    fmpq_mul_fmpz(st_data->f+i, st_data->f+i, st_data->modlist+i);
  }

  fmpz_mat_init(st_data->binom_mat, d+1, d+1);
  for (i=0; i<=d; i++)
    for (j=0; j<=d; j++)
      fmpz_bin_uiui(fmpq_mat_entry(st_data->binom_mat, i, j), i, j);
  
  st_data->sum_mats = (fmpq_mat_t *)malloc((d+1)*sizeof(fmpq_mat_t));
  for (i=0; i<=d; i++) {

    fmpq_mat_init(st_data->sum_mats[i], 9, d+1);
    fmpq_mat_zero(st_data->sum_mats[i]);

    arith_chebyshev_t_polynomial(pol, i);
    for (j=0; j<=d; j++) {
      
      /* Row 0: coeffs of 2*(i-th Chebyshev polynomial)(x/2). */
      if (j <= i) {
	k1 = fmpq_mat_entry(st_data->sum_mats[i], 0, j);
	fmpq_set_fmpz_frac(k1, fmpz_poly_get_coeff_ptr(pol, j), const1);
	fmpz_mul_2exp(m, const1, j);
	fmpq_div_fmpz(k1, k1, m);
	fmpz_set_ui(m, 2);
	fmpq_mul_fmpz(k1, k1, m);
      }

      /* Row 1: coeffs of row 0 from matrix i-2, multiplied by -2. */
      if (i >= 2) {
	k1 = fmpq_mat_entry(st_data->sum_mats[i], 1, j);
	fmpq_set(k1, fmpq_mat_entry(st_data->sum_mats[i-2], 0, j));
	fmpz_set_si(m, -2);
	fmpq_mul_fmpz(k1, k1, m);
      }

      /* Row 2: coeffs of row 0 from matrix i-2, shifted by 2. */
      if (i>= 2 && j >= 2) {
	k1 = fmpq_mat_entry(st_data->sum_mats[i], 2, j);
	fmpq_set(k1, fmpq_mat_entry(st_data->sum_mats[i-2], 0, j-2));
      }

      /* Row 3: coeffs of (2+x)^i. */
      if (j<= i) {
	k1 = fmpq_mat_entry(st_data->sum_mats[i], 3, j);
	fmpq_set_fmpz_frac(k1, fmpq_mat_entry(st_data->binom_mat, i, j), const1);
	fmpq_mul_2exp(k1, k1, i-j);
      }
      
      /* Row 4: coeffs of (2+x)^(i-1). */
      if (i >= 1) {
	k1 = fmpq_mat_entry(st_data->sum_mats[i], 4, j);
	fmpq_set(k1, fmpq_mat_entry(st_data->sum_mats[i-1], 3, j));	
      }

      /* Row 5: coeffs of (2+x)^(i-2). */
      if (i>=2)	{
	k1 = fmpq_mat_entry(st_data->sum_mats[i], 5, j);
	fmpq_set(k1, fmpq_mat_entry(st_data->sum_mats[i-2], 3, j));	
      }

      /* Row 6: coeffs of (-2+x)^i. */
      k1 = fmpq_mat_entry(st_data->sum_mats[i], 6, j);
      fmpq_set(k1, fmpq_mat_entry(st_data->sum_mats[i], 3, j));
      if ((i-j)%2==1) fmpq_neg(k1, k1);

      /* Row 7: coeffs of (-2+x)^(i-1). */
      if (i >= 1) {
	k1 = fmpq_mat_entry(st_data->sum_mats[i], 7, j);
	fmpq_set(k1, fmpq_mat_entry(st_data->sum_mats[i-1], 6, j));	
      }

      /* Row 8: coeffs of (-2+x)^(i-2). */
      if (i >= 2) {
	k1 = fmpq_mat_entry(st_data->sum_mats[i], 8, j);
	fmpq_set(k1, fmpq_mat_entry(st_data->sum_mats[i-2], 6, j));
      }
    }
  }
  
  fmpz_poly_clear(pol);
  fmpz_clear(m);
  fmpz_clear(const1);

  return(st_data);
}

ps_dynamic_data_t *ps_dynamic_init(int d, int *Q0) {
  ps_dynamic_data_t *dy_data;
  int i;

  dy_data = (ps_dynamic_data_t *)malloc(sizeof(ps_dynamic_data_t));
  dy_data->d = d;

  /* Initialize mutable quantities */
  dy_data->n = d;
  dy_data->count = 0;
  dy_data->ascend = 0;
  dy_data->pol = _fmpz_vec_init(d+1);
  if (Q0 != NULL) 
    for (i=0; i<=d; i++) 
      fmpz_set_si(dy_data->pol+i, Q0[i]);
  
  fmpq_mat_init(dy_data->sum_col, d+1, 1);
  fmpq_set_si(fmpq_mat_entry(dy_data->sum_col, 0, 0), d, 1);

  dy_data->upper = _fmpz_vec_init(d+1);

  /* Allocate scratch space */
  fmpq_mat_init(dy_data->sum_prod, 9, 1);
  dy_data->wlen = 4*d+12;
  dy_data->w = _fmpz_vec_init(dy_data->wlen);
  dy_data->w2len = 5;
  dy_data->w2 = _fmpq_vec_init(dy_data->w2len);
  return(dy_data);
}

ps_dynamic_data_t *ps_dynamic_clone(ps_dynamic_data_t *dy_data) {
  ps_dynamic_data_t *dy_data2;
  int i, d = dy_data->d;

  dy_data2 = ps_dynamic_init(d, NULL);
  dy_data2->n = dy_data->n;
  dy_data2->count = dy_data->count;
  dy_data2->ascend = dy_data->ascend;
  _fmpz_vec_set(dy_data2->pol, dy_data->pol, d+1);
  _fmpz_vec_set(dy_data2->upper, dy_data->upper, d+1);
  fmpq_mat_set(dy_data2->sum_col, dy_data->sum_col);
  return(dy_data2);
}

ps_dynamic_data_t *ps_dynamic_split(ps_dynamic_data_t *dy_data) {
  if (dy_data==NULL) return(NULL);

  ps_dynamic_data_t *dy_data2;
  int i, d = dy_data->d, n = dy_data->n;

  for (i=d; i>n+1; i--)
    if (fmpz_cmp(dy_data->pol+i, dy_data->upper+i) <0) {
      
      dy_data2 = ps_dynamic_clone(dy_data);
      fmpz_set(dy_data->upper+i, dy_data->pol+i);
      dy_data2->n = i-1;
      dy_data2->ascend = 1;
      return(dy_data2);
  }
  return(NULL);
}

void extract_pol(int *Q, ps_dynamic_data_t *dy_data) {
  int i;
  fmpz *pol = dy_data->pol;
  for (i=0; i<=dy_data->d; i++)
    Q[i] = fmpz_get_si(pol+i);
}

int extract_count(ps_dynamic_data_t *dy_data) {
  return(dy_data->count);
}

void ps_static_clear(ps_static_data_t *st_data) {
  int i, d = st_data->d;
  fmpz_clear(st_data->a);
  fmpz_clear(st_data->b);
  fmpz_mat_clear(st_data->binom_mat);
  _fmpq_vec_clear(st_data->f, d+1);
  _fmpz_vec_clear(st_data->modlist, d+1);
  for (i=0; i<=d; i++) 
    fmpq_mat_clear(st_data->sum_mats[i]);
  free(st_data->sum_mats);
  free(st_data);
}

void ps_dynamic_clear(ps_dynamic_data_t *dy_data) {
  int d = dy_data->d;
  _fmpz_vec_clear(dy_data->pol, d+1);
  _fmpz_vec_clear(dy_data->upper, d+1);
  fmpq_mat_clear(dy_data->sum_col);
  fmpq_mat_clear(dy_data->sum_prod);
  _fmpz_vec_clear(dy_data->w, dy_data->wlen);
  _fmpq_vec_clear(dy_data->w2, dy_data->w2len);
  free(dy_data);
}

/* Return values: 
   1: if lower <= upper
   0: otherwise. 
   Both cases include the option n=0, in which case we simply check
   admissibility of the polynomial (there being no further coefficients
   to control).
*/
int set_range_from_power_sums(ps_static_data_t *st_data,
			  ps_dynamic_data_t *dy_data) {
  int i, j, r, r1, r2;
  int d = st_data->d;
  int n = dy_data->n;
  int k = d+1-n;
  fmpz *modulus = st_data->modlist + n-1; // Cannot dereference until n>0
  fmpz *pol = dy_data->pol;
  fmpq *f;
    
  /* Allocate temporary variables from persistent scratch space. */
  fmpz *tpol = dy_data->w;

  fmpz *t0z = dy_data->w + d + 1;
  fmpz *lower = dy_data->w + d + 2;
  fmpz *upper = dy_data->w + d + 3;
  
  fmpq *t0q = dy_data->w2;
  fmpq *t1q = dy_data->w2 + 1;
  fmpq *t2q = dy_data->w2 + 2;
  fmpq *t3q = dy_data->w2 + 3;

  /* Subroutines to adjust lower and upper bounds. */

  void set_lower(const fmpq_t val) {
    fmpq_div(t0q, val, f);
    fmpq_ceil(lower, t0q);
  }
  
  void set_upper(const fmpq_t val) {
    fmpq_div(t0q, val, f);
    fmpq_floor(upper, t0q);
  }

  void change_lower(const fmpq_t val) {
    fmpq_div(t0q, val, f);
    fmpq_ceil(t0z, t0q);
    if (fmpz_cmp(t0z, lower) > 0) fmpz_set(lower, t0z);
  }
  
  void change_upper(const fmpq_t val) {
    fmpq_div(t0q, val, f);
    fmpq_floor(t0z, t0q);
    if (fmpz_cmp(t0z, upper) < 0) fmpz_set(upper, t0z);
  }
    
  /* Compute the divided n-th derivative of pol, 
     then determine whether its roots are all in [-2, 2]. */
  for (i=0; i<=k-1; i++)
    fmpz_mul(tpol+i, fmpz_mat_entry(st_data->binom_mat, n+i, n), pol+n+i);

  r = _fmpz_poly_all_roots_in_interval(tpol, k, st_data->a, st_data->b, 
				       dy_data->w+d+1);
  /* If r=0, abort. 
     If r=1 and k>d, return trivially; no further coefficients to find.
     If r=1 and k<=d, continue to compute bounds.
  */

  if (r==0) return(0);
  if ((k > d) || fmpz_is_zero(modulus)) return(1);

  /* Compute the k-th power sum. */
  f = fmpq_mat_entry(dy_data->sum_col, k, 0);
  fmpq_set_si(f, -k, 1);
  fmpq_mul_fmpz(f, f, pol+d-k);
  fmpq_div_fmpz(f, f, pol+d);
  for (i=1; i<k; i++) {
    fmpq_set_fmpz_frac(t0q, pol+d-i, pol+d);
    fmpq_neg(t0q, t0q);
    fmpq_addmul(f, t0q, fmpq_mat_entry(dy_data->sum_col, k-i, 0));
  }
  
  /* Initialize bounds using power sums of the asymmetrized polynomial. */
  f = st_data->f + n-1;
  fmpq_mat_mul(dy_data->sum_prod, st_data->sum_mats[k], dy_data->sum_col);

  fmpq_set_si(t1q, 2*d, 1);
  fmpq_sub(t0q, fmpq_mat_entry(dy_data->sum_prod, 0, 0), t1q);
  set_lower(t0q);
  fmpq_add(t0q, fmpq_mat_entry(dy_data->sum_prod, 0, 0), t1q);
  set_upper(t0q);

  /* Apply Descartes' rule of signs at -2, +2. */
  /* Currently tpol is the divided n-th derivative of pol;
     evaluate at the endpoints. */
  _fmpz_poly_evaluate_fmpz(t0z, tpol, k, st_data->a);
  r1 = fmpz_sgn(t0z);
  _fmpz_poly_evaluate_fmpz(t0z, tpol, k, st_data->b);
  r2 = fmpz_sgn(t0z);
  
  /* Now set tpol to the divided (n-1)-st derivative of pol.
     then evaluate again. */
  for (i=k; i>=1; i--) {
    fmpz_mul_si(tpol+i, tpol+i-1, n);
    fmpz_divexact_si(tpol+i, tpol+i, i);
  }
  fmpq_set_si(t2q, -k, 1);
  fmpq_div_fmpz(t2q, t2q, pol+d);
  
  fmpz_set(tpol, pol+d-k);    
  _fmpz_poly_evaluate_fmpz(t0z, tpol, k+1, st_data->a);
  fmpq_mul_fmpz(t1q, t2q, t0z);
  if (r1 >= 0) change_upper(t1q);
  if (r1 <= 0) change_lower(t1q);
  
  _fmpz_poly_evaluate_fmpz(t0z, tpol, k+1, st_data->b);
  fmpq_mul_fmpz(t1q, t2q, t0z);
  if (r2 >= 0) change_lower(t1q);
  if (r2 <= 0) change_upper(t1q);
  
  /* Additional bounds based on power sums. */
  if ((fmpz_cmp(lower, upper) <= 0) && k >= 2) {
    fmpq_add(t1q, fmpq_mat_entry(dy_data->sum_prod, 1, 0),
	     fmpq_mat_entry(dy_data->sum_prod, 2, 0));
    fmpq_set_si(t2q, 4*d, 1);
    fmpq_sub(t0q, t1q, t2q);
    // fmpq_sub_si(t3q, t1q, 4*d);
    if (k==2) fmpq_div_fmpz(t0q, t0q, st_data->b); // b=2
    change_lower(t0q);
    fmpq_add(t0q, t1q, t2q);
    // fmpq_add_si(t3q, t1q, 4*d);
    if (k==2) fmpq_div_fmpz(t0q, t0q, st_data->b); // b=2
    change_upper(t0q);
    
    /* t1q, t2q, t3q are no longer needed, so can be reassigned. */
    t1q = fmpq_mat_entry(dy_data->sum_prod, 3, 0);
    t2q = fmpq_mat_entry(dy_data->sum_prod, 4, 0);
    t3q = fmpq_mat_entry(dy_data->sum_prod, 5, 0);
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
    
    t3q = fmpq_mat_entry(dy_data->sum_prod, 8, 0);
    if ((k%2 == 0) && (fmpq_sgn(t3q) > 0)) {
      t2q = fmpq_mat_entry(dy_data->sum_prod, 7, 0);
      fmpq_mul(t0q, t2q, t2q);
      t1q = fmpq_mat_entry(dy_data->sum_prod, 6, 0);
      fmpq_div(t0q, t0q, t3q);
      fmpq_sub(t0q, t1q, t0q);
      change_upper(t0q);
    } else if ((k%2 == 1) && (fmpq_sgn(t3q) < 0)) {
      t2q = fmpq_mat_entry(dy_data->sum_prod, 7, 0);
      fmpq_mul(t0q, t2q, t2q);
      t1q = fmpq_mat_entry(dy_data->sum_prod, 6, 0);
      fmpq_div(t0q, t0q, t3q);
      fmpq_sub(t0q, t1q, t0q);
      change_lower(t0q);
    }
    fmpq_set_si(t0q, 4, 1);
    fmpq_mul(t0q, t0q, t2q);
    // fmpq_mul_si(t0q, t2q, 4);
    fmpq_add(t0q, t0q, t1q);
    if (k%2 == 0) change_lower(t0q);
    else change_upper(t0q);
    
    if (k%2 == 0) {
	fmpq_set_si(t0q, -4, 1);
	fmpq_mul(t0q, t0q, fmpq_mat_entry(dy_data->sum_col, k-2, 0));
      // fmpq_mul_si(t0q, fmpq_mat_entry(dy_data->sum_col, k-2, 0), -4);
	fmpq_add(t0q, t0q, fmpq_mat_entry(dy_data->sum_col, k, 0));
	change_lower(t0q);	
      }
  }

  if (fmpz_cmp(lower, upper) > 0) return(0);

  /* Set the new upper bound. */
  fmpz_mul(upper, upper, modulus);
  fmpz_add(dy_data->upper+n-1, pol+n-1, upper);
  /* Correct the k-th power sum. */
  t1q = fmpq_mat_entry(dy_data->sum_col, k, 0);
  fmpq_mul_fmpz(t0q, f, lower);
  fmpq_sub(t1q, t1q, t0q);
  /* Set the new polynomial value. */
  fmpz_mul(lower, lower, modulus);
  fmpz_add(pol+n-1, pol+n-1, lower);

  return(1);

}

/* Return values:
    1: if a solution has been found
    0: if the tree has been exhausted
   -1: if the maximum number of nodes has been reached
*/

int next_pol(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data) {
  if (dy_data==NULL) return(0);

  int d = st_data->d;
  int verbosity = st_data->verbosity;
  int node_count = st_data->node_count;
  fmpz *modlist = st_data->modlist;

  int ascend = dy_data->ascend;
  int n = dy_data->n;
  int count = dy_data->count;
  fmpz *upper = dy_data->upper;
  fmpz *pol = dy_data->pol;

  int i, t;
  fmpq *tq;

  if (n>d) return(0);
  while (1) {
    if (ascend) {
      n += 1;
      if (n>d) { t=0; break; }
    } else {
      if (d-n <= verbosity) {
	_fmpz_vec_print(pol+n, d-n+1);
	printf("\n");
      }
      count += 1;
      if (node_count != -1 && count >= node_count) { t= -1; break; }
      dy_data->n = n;
      if (set_range_from_power_sums(st_data, dy_data)) {
	  n -= 1;
	  if (n<0) { t=1; break; }
	  continue;
	}
    }
    if (fmpq_is_zero(modlist+n)) ascend = 1;
    else {
      fmpz_add(pol+n, pol+n, modlist+n);
      if (fmpz_cmp(pol+n, upper+n) > 0) ascend = 1;
      else {
	ascend = 0;
	/* Update the (d-n)-th power sum. */
	tq = fmpq_mat_entry(dy_data->sum_col, d-n, 0);
	fmpq_sub(tq, tq, st_data->f+n);
      }
    }
  }
  dy_data->ascend = (n<0);
  dy_data->n = n;
  dy_data->count = count;
  return(t);
}
