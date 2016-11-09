#ifndef KSKIPBICG_H_INCLUDED__
#define KSKIPBICG_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include "../share.h"
#include "../io.h"
#include "../tools.h"
#include "../blas.h"


void KSKIPBICG_init(double **Ap, double **Ar, double *theta, double *eta, double *rho, double *phi, double *rvec, double *pvec, double *r_vec, double *p_vec, double *Av, double *xvec, const int N, const int i_kskip);

int KSKIPBICG_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, const struct Parameter *para, const int N, const int NNZ, const bool f_isinner);


#endif //KSKIPBICG_H_INCLUDED__

