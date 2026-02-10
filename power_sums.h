#pragma once
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_poly.h>

#if defined(_OPENMP)
  #include <omp.h>
#endif

typedef struct ps_static_data {
  int d, force_squarefree;
  long node_limit;
  fmpz_t q;
  fmpz *modlist, *binom_mat, *sum_mats, *eval_pm2_mats;
} ps_static_data_t;

typedef struct ps_dynamic_data {
  int d, n, ascend, flag;
  long node_count;
  fmpz *pol, *sympol, *upper, *power_sums_num, *hankel_dets;

  /* Scratch space */
  fmpz *w;
  long wlen; /* = 3*d+10 */
} ps_dynamic_data_t;

int has_openmp();
int num_threads();
int is_mpz(fmpz f);

ps_static_data_t *ps_static_init(int d, fmpz_t q, fmpz_t lead, fmpz *modlist, 
				 long node_limit, int force_squarefree);
ps_dynamic_data_t *ps_dynamic_init(int d, fmpz_t q, fmpz *coefflist);
void ps_static_clear(ps_static_data_t *st_data);
void ps_dynamic_clear(ps_dynamic_data_t *dy_data);

void ps_dynamic_split(ps_dynamic_data_t *dy_data, ps_dynamic_data_t *dy_data2);
void next_pol(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int max_steps);
