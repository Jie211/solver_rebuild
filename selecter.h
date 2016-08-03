#ifndef SELECTER_H_INCLUDED__
#define SELECTER_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "io.h"
#include "share.h"
#include "./CRS/cg.h"
#include "./CRS/cr.h"
#include "./CRS/gcr.h"
#include "./CRS/gmres.h"
#include "./CRS/kskipcg.h"
#include "./CRS/kskipcr.h"

int outer_selecter(struct Parameter *para, double *bvec, double *xvec, double *val, int *col, int *ptr, const int N, const int NNZ);

#endif //SELECTER_H_INCLUDED__

