#ifndef ALL_ROOTS_IN_INTERVAL
#define ALL_ROOTS_IN_INTERVAL

#include <fmpz_poly.h>

int _fmpz_poly_all_roots_in_interval(fmpz *poly, slong n, 
                                     fmpz const * a, fmpz const * b, fmpz *w);
int _fmpz_poly_all_roots_real(fmpz *poly, slong n, fmpz *w);

#endif

