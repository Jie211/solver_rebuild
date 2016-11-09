#ifndef CR_H_INCLUDED__
#define CR_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include "../share.h"
#include "../io.h"
#include "../tools.h"
#include "../blas.h"

void CR_init(double *rvec, double *pvec, double *qvec, double *svec, double *xvec, const int N);

int CR_CRS(double *val, int *col, int *ptr, double *Tval, int *Tcol, int *Tptr, double *bvec, double *xvec, const struct Parameter *para, const int N, const int NNZ, const bool f_isinner);

#endif //CR_H_INCLUDED__

