#ifndef START_H_INCLUDED__
#define START_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "io.h"
#include "share.h"
#include "tools.h"
#include "selecter.h"
#include "blas.h"
#include "cuda_func.cuh"

int csr_start(int argc, char *argv[]);

#endif //START_H_INCLUDED__

