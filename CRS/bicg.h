#ifndef BICG_H_INCLUDED__
#define BICG_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include "../share.h"
#include "../io.h"
#include "../tools.h"
#include "../blas.h"

void BICG_Init(double *rvec, double *r_vec, double *pvec, double *p_vec, double *Av, double *xvec, const int N);

int BICG_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, const struct Parameter *para, const int N, const int NNZ, const bool f_isinner);

#endif //BICG_H_INCLUDED__

