#ifndef IO_H_INCLUDED__
#define IO_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdbool.h>

#include "share.h"
#include "tools.h"

void init_ver(struct Parameter *para);
int get_opt(int argc, char *argv[], struct Parameter *para);
int check_opt(struct Parameter *para);
int check_solver(char *optarg, enum SolverName *solver);
void show_opt(struct Parameter *para);

#endif //IO_H_INCLUDED__

