#ifndef VPCG_H_INCLUDED__
#define VPCG_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include "../share.h"
#include "../io.h"
#include "../tools.h"
#include "../blas.h"
#include "../selecter.h"

void VPCG_Init(double *rvec, double *pvec, double *Av, double *zvec, double *xvec, const int N);

int VPCG_CRS(double *val, int *col, int *ptr, double *Tval, int *Tcol, int *Tptr, double *bvec, double *xvec, struct Parameter *para, const int N, const int NNZ, const bool f_isinner);
#endif //VPCG_H_INCLUDED__

