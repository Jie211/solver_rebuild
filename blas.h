#ifndef BLAS_H_INCLUDED__
#define BLAS_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double norm_1_d(double *v, const int N);

double norm_2_d(double *v, const int N);

void MV_mult_CSR(double *out, const double *val, const int *col, const int *ptr, const double *vec, const int N);

void vec_sub(double *out, const double *x, const double *y, const int N);

void vec_add(double *out, const double *x, const double *y, const int N);

double dot_d(const double *x, const double *y, const int N);

void scalar_d(double *out, const double a, const double *x, const int N);

void scalar_xpy_d(double *out, const double a, const double *x, const double *y, const int N);

double error_d_CRS(double *val, const int *col, const int *ptr, const double *b, const double *x_new, const double *x_0, const int N);

void solve_Hye(double *h, double *y, double *e, const int n, const int size);

void cal_arap_kskipcg_d(double **Ar, double **Ap, double *val, int *col, int *ptr, double *rvec, double *pvec, const int N, const int kskip);

void cal_deltaetazeta_kskipcg_d(double *delta, double *eta, double *zeta, double **Ar, double **Ap, double *rvec, double *pvec, const int N, const int kskip);

#endif //BLAS_H_INCLUDED__

