#include <fmpz_poly.h>
#include <fmpq.h>
#include <fmpq_mat.h>

typedef struct power_sums_data {
  int d, lead;
  fmpq_mat_t *sum_mats;
  fmpq_mat_t sum_col, sum_prod;
} power_sums_data_t;

power_sums_data_t *ranger_init(int d, int lead, int *modlist);
void ranger_clear(power_sums_data_t *data);
int range_from_power_sums(int *bounds, power_sums_data_t *data,
			  int *pol, int k);

