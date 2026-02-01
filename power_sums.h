#pragma once
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>

typedef struct ps_static_data {
  int d, sign, force_squarefree, q_is_1;
  long node_limit;
  fmpz_t a, b, lead, q;
  fmpz *modlist, *binom_mat, *sum_mats, *eval_pm2_mats;
  fmpq *f;
} ps_static_data_t;

typedef struct ps_dynamic_data {
  int d, n, ascend, flag;
  long node_count;
  fmpq_mat_t power_sums, hankel_mat, hankel_dets[2];
  fmpz *pol, *sympol, *upper;

  /* Scratch space */
  fmpz *w;
  long wlen; /* = 3*d+8 */
  fmpq *w2;
  long w2len; /* = 4 */
} ps_dynamic_data_t;

int has_openmp();
int num_threads();
ps_static_data_t *ps_static_init(int d, fmpz_t q, fmpz_t lead, fmpz *modlist, 
				 long node_limit, int force_squarefree);
ps_dynamic_data_t *ps_dynamic_init(int d, fmpz_t q, fmpz *coefflist);
void ps_static_clear(ps_static_data_t *st_data);
void ps_dynamic_clear(ps_dynamic_data_t *dy_data);
void extract_pol(int *Q, ps_dynamic_data_t *dy_data);
void ps_dynamic_split(ps_dynamic_data_t *dy_data, ps_dynamic_data_t *dy_data2);
inline void step_forward(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int n, 
                         fmpz_t step);
void next_pol(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int max_steps);
