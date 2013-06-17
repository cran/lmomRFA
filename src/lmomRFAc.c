#include <R.h>

void F77_SUB(ranget)(void) { GetRNGstate(); }

void F77_SUB(ranput)(void) { PutRNGstate(); }

void F77_SUB(curand)(int *n, double *x) {
  int i;
  for (i = 0; i < *n; i++)
	   x[i] = unif_rand();
}
