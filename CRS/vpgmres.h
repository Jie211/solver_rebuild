#ifndef VPGMRES_H_INCLUDED__
#define VPGMRES_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include "../share.h"
#include "../io.h"
#include "../tools.h"
#include "../blas.h"
#include "../selecter.h"

void VPGMRES_Init(double *rvec, double *axvec, double *evec, double *vvec, double *vmtx, double *hmtx, double *yvec, double *wvec, double *avvec, double *hvvec, double *cvec, double *svec, double *x0vec, double *tmpvec, double *zmtx, double *zvec, double *xvec, const int N, const int i_restart);

int VPGMRES_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, struct Parameter *para, const int N, const int NNZ, const bool f_isinner);

#endif //VPGMRES_H_INCLUDED__

