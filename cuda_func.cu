#include "cuda_func.h"

void malloc_cuda_1d(int num_x, double *arr)
{
  cudaMalloc((void**)arr, num_x*sizeof(double));
}
