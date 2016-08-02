#ifndef TOOLS_H_INCLUDED__
#define TOOLS_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <omp.h>

#include "share.h"

void error_log(char *output);

void warning_log(char *output);

void normal_log(char *output);

double *malloc_1d(int num);

double **malloc_2d(int num_x, int num_y);

int *malloc_1i(int num);

int** malloc_2i(int num_x, int num_y);

void free_1d(double *ptr);

void free_2d(double **ptr, int num);

void free_1i(int *ptr);

void free_2i(int **ptr, int num);

int error_handle(int error_code, char *msg);

int set_openmp_thread(const int thread);

FILE *file_init(const char *name, const char *mode);

void vec_init(double *v, const double val, const int size);

void vec_init_2(double *v, const double val, const int sizex, const int sizey);

void vec_copy(double *t, double *f, const int size);

void file_print(FILE *fp, double *vec, const int N);

double error_d_CRS(double *val, const int *col, const int *ptr, const double *b, const double *x_new, const double *x_0, const int N);
void vec_add(double *out, double *x, double *y, const int N);

#endif //TOOLS_H_INCLUDED__

