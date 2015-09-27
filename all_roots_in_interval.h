#ifndef ALL_ROOTS_IN_INTERVAL
#define ALL_ROOTS_IN_INTERVAL

#include <fmpz_poly.h>

int _fmpz_poly_all_roots_in_interval(fmpz *poly, slong n, 
                                     fmpz const * a, fmpz const * b, fmpz *w);
int fmpz_poly_all_roots_in_fixed_interval(int *pol, int deg, int k);

#endif

