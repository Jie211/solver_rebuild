#ifndef GCR_H_INCLUDED__
#define GCR_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include "../share.h"
#include "../io.h"
#include "../tools.h"
#include "../blas.h"


void GCR_init(double *rvec, double *Av, double *x_0, double *qq, double **qvec, double **pvec, const int N, const int restart);


int GCR_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, const struct Parameter *para, const int N, const int NNZ, const bool f_isinner);


#endif //GCR_H_INCLUDED__
