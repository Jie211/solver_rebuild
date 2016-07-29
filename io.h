#ifndef IO_H_INCLUDED__
#define IO_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <getopt.h>
#include <stdbool.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "share.h"
#include "tools.h"

void init_ver(struct Parameter *para);
int get_opt(int argc, char *argv[], struct Parameter *para);
int check_opt(struct Parameter *para);
int check_solver(char *optarg, enum SolverName *solver);
void show_opt(struct Parameter *para);
int find_mat(struct Parameter *para);
int get_mat_head(const struct Parameter *para, int *N, int *NNZ);
int get_mat_data(const struct Parameter *para, int *col, int *ptr, double *val, double *bvec, double *xvec, const int N, const int NNZ);
#endif //IO_H_INCLUDED__

