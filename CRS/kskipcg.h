#ifndef KSKIPCG_H_INCLUDED__
#define KSKIPCG_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include "../share.h"
#include "../io.h"
#include "../tools.h"
#include "../blas.h"

void KSKIPCG_init(double **Ap, double **Ar, double *delta, double *eta, double *zeta, double *rvec, double *pvec, double *Av, double *xvec, const int N, const int kskip);
int KSKIPCG_CRS(double *val, int *col, int *ptr, double *Tval, int *Tcol, int *Tptr, double *bvec, double *xvec, const struct Parameter *para, const int N, const int NNZ, const bool f_isinner);
#endif //KSKIPCG_H_INCLUDED__

