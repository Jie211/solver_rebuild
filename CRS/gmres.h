#ifndef GMRES_H_INCLUDED__
#define GMRES_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include "../share.h"
#include "../io.h"
#include "../tools.h"
#include "../blas.h"

int GMRES_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, const struct Parameter *para, const int N, const int NNZ, const bool isinner);

void GMRES_init(double *rvec, double *axvec, double *evec, double *vvec, double *vmtx, double *hmtx, double *yvec, double *wvec, double *avvec, double *hvvec, double *cvec, double *svec, double *x0vec, double *tmpvec, double *x_0, double *xvec, const int N, const int restart);
#endif //GMRES_H_INCLUDED__

