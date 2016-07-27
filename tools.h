#ifndef TOOLS_H_INCLUDED__
#define TOOLS_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

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

void init_ver(char *c_matrix, 
    char *c_outer_solver, int *i_outer_maxloop, double *d_outer_eps, int *i_outer_restart, int *i_outer_kskip, int *i_outer_fix,
    char *c_inner_solver, int *i_inner_maxloop, double *d_inner_eps, int *i_inner_restart, int *i_inner_kskip, int *i_inner_fix,
   int *i_thread, bool *f_cuda, bool *f_verbose);

#endif //TOOLS_H_INCLUDED__

