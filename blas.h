#ifndef BLAS_H_INCLUDED__
#define BLAS_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double norm_1_d(double *v, const int N);

double norm_2_d(double *v, const int N);

void MV_mult_CSR(double *out, const double *val, const int *col, const int *ptr, const double *vec, const int N);

void vec_sub(double *out, const double *x, const double *y, const int N);
double dot_d(const double *x, const double *y, const int N);

void scalar_d(double *out, const double a, const double *x, const int N);
void scalar_xpy_d(double *out, const double a, const double *x, const double *y, const int N);

double error_d_CRS(double *val, const int *col, const int *ptr, const double *b, const double *x_new, const double *x_0, const int N);

#endif //BLAS_H_INCLUDED__

