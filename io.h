#ifndef IO_H_INCLUDED__
#define IO_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdbool.h>

#include "tools.h"

void init_ver(char *c_matrix, 
    char *c_outer_solver, int *i_outer_maxloop, double *d_outer_eps, int *i_outer_restart, int *i_outer_kskip, int *i_outer_fix,
    char *c_inner_solver, int *i_inner_maxloop, double *d_inner_eps, int *i_inner_restart, int *i_inner_kskip, int *i_inner_fix,
   int *i_thread, bool *f_cuda, bool *f_verbose);

int get_cmd(int argc, char *argv[], 
    char *c_matrix, 
    char *c_outer_solver, int *i_outer_maxloop, double *d_outer_eps, int *i_outer_restart, int *i_outer_kskip, int *i_outer_fix,
    char *c_inner_solver, int *i_inner_maxloop, double *d_inner_eps, int *i_inner_restart, int *i_inner_kskip, int *i_inner_fix,
   int *i_thread, bool *f_cuda, bool *f_verbose);

int check_cmd(char *c_matrix, 
    char *c_outer_solver, int *i_outer_maxloop, double *d_outer_eps, int *i_outer_restart, int *i_outer_kskip, int *i_outer_fix,
    char *c_inner_solver, int *i_inner_maxloop, double *d_inner_eps, int *i_inner_restart, int *i_inner_kskip, int *i_inner_fix,
   int *i_thread, bool *f_cuda, bool *f_verbose);

void show_cmd(char *c_matrix, 
    char *c_outer_solver, int *i_outer_maxloop, double *d_outer_eps, int *i_outer_restart, int *i_outer_kskip, int *i_outer_fix,
    char *c_inner_solver, int *i_inner_maxloop, double *d_inner_eps, int *i_inner_restart, int *i_inner_kskip, int *i_inner_fix,
   int *i_thread, bool *f_cuda, bool *f_verbose);

#endif //IO_H_INCLUDED__

