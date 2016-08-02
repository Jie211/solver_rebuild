#ifndef CG_H_INCLUDED__
#define CG_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include "../share.h"
#include "../io.h"
#include "../tools.h"
#include "../blas.h"

void CG_Init(double *rvec, double *pvec, double *Av, double *xvec, const int N);

int CG_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, const struct Parameter *para, const int N, const int NNZ, const bool f_isinner);
#endif //CG_H_INCLUDED__

