#ifndef VPGCR_H_INCLUDED__
#define VPGCR_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include "../share.h"
#include "../io.h"
#include "../tools.h"
#include "../blas.h"
#include "../selecter.h"

void VPGCR_Init(double *rvec, double *zvec, double *Av, double **qvec, double **pvec, double *qq, double *xvec, const int N, const int restart);

int VPGCR_CRS(double *val, int *col, int *ptr, double *Tval, int *Tcol, int *Tptr, double *bvec, double *xvec, struct Parameter *para, const int N, const int NNZ, const bool f_isinner);
#endif //VPGCR_H_INCLUDED__

