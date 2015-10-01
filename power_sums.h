#include <fmpz_poly.h>
#include <fmpq.h>
#include <fmpq_mat.h>

typedef struct ps_static_data {
  int d, lead, verbosity, node_count;
  fmpz_t a, b;
  fmpz_mat_t binom_mat;
  fmpz *modlist;
  fmpq_mat_t *sum_mats;
  fmpq_t *f;
} ps_static_data_t;

typedef struct ps_dynamic_data {
  int d, n, count;
  fmpq_mat_t sum_col, sum_prod;
  fmpz *pol, *upper;

  /* Scratch space */
  fmpz *w;
  slong wlen; /* = 4*d+12 */
  fmpq *w2;
  slong w2len; /* = 5 */
} ps_dynamic_data_t;

ps_static_data_t *ps_static_init(int d, int lead, int *modlist, 
				 int verbosity, int node_count);
ps_dynamic_data_t *ps_dynamic_init(int d, int *Q0);
void ps_static_clear(ps_static_data_t *st_data);
void ps_dynamic_clear(ps_dynamic_data_t *dy_data);
int next_pol(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int *Q);

