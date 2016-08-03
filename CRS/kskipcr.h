#ifndef KSKIPCR_H_INCLUDED__
#define KSKIPCR_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include "../share.h"
#include "../io.h"
#include "../tools.h"
#include "../blas.h"

void KSKIPCR_init(double *Ar, double *Ap, double *delta, double *eta, double *zeta, double *rvec, double *pvec, double *Av, double *xvec, const int N, const int kskip);

int KSKIPCR_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, const struct Parameter *para, const int N, const int NNZ, const bool f_isinner);

#endif //KSKIPCR_H_INCLUDED__

