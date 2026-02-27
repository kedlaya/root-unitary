#include <stdio.h>
#include <stdlib.h>
#include "../power_sums.c"

int main(int argc, char* argv[]) {
  int i, d, d0, flag = 1;
  long q, lead, count = 0;
  fmpz_t temp_lead, temp_q;
  fmpz *temp_array;
  ps_static_data_t *st_data;
  ps_dynamic_data_t *dy_data;
  
  if (argc != 4) {
    fprintf(stderr, "Format: weilpoly d q lead\n");
    return(1);
  }
  d0 = atoi(argv[1]);
  q = atoi(argv[2]);
  lead = atoi(argv[3]);
  fprintf(stderr, "Computing Weil polynomials with d = %d, q = %ld, lead = %ld\n", d0, q, lead);

  d = d0/2;
  fmpz_init(temp_lead);
  fmpz_set_si(temp_lead, lead);
  fmpz_init(temp_q);
  fmpz_set_si(temp_q, q);
  temp_array = _fmpz_vec_init(d+1);
  fmpz_zero(temp_array+d);
  for (i=1; i<=d; i++) fmpz_one(temp_array+i);
  st_data = ps_static_init(d, temp_q, temp_lead, temp_array, -1, 0);
  _fmpz_vec_zero(temp_array, d);
  fmpz_set(temp_array+d, temp_lead);
  dy_data = ps_dynamic_init(d, temp_q, temp_array);
  
  while (flag) {
    flag = next_pol(st_data, dy_data, 1000);
    if (flag == 2) {
      for (i=0; i<=d0; i++) { 
        fmpz_print(dy_data->sympol+i);
        if (i<d0) printf(" ");
      }
      printf("\n");
      count++;
    }
  }
  
  ps_static_clear(st_data);
  ps_dynamic_clear(dy_data);
  fprintf(stderr, "Total number of polynomials: %ld\n", count);
  return(0);
}
