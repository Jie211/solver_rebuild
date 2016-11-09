#ifndef VPCR_H_INCLUDED__
#define VPCR_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include "../share.h"
#include "../io.h"
#include "../tools.h"
#include "../blas.h"
#include "../selecter.h"

void VPCR_Init(double *rvec, double *pvec, double *zvec, double *Av, double *Ap, double *xvec, const int N);

int VPCR_CRS(double *val, int *col, int *ptr, double *Tval, int *Tcol, int *Tptr, double *bvec, double *xvec, struct Parameter *para, const int N, const int NNZ, const bool f_isinner);

#endif //VPCR_H_INCLUDED__

