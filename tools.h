#ifndef TOOLS_H_INCLUDED__
#define TOOLS_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>

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

void error_handle(int error_code, char *msg);

#endif //TOOLS_H_INCLUDED__

