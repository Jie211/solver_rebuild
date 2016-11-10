#include "cuda_func.cuh"

void malloc_cuda_1d(int num_x, double *arr)
{
  cudaMalloc((void**)&arr, num_x*sizeof(double));
}

void malloc_cuda_1i(int num_x, int *arr)
{
  cudaMalloc((void**)&arr, num_x*sizeof(int));
}

void free_cuda_1d(double *arr)
{
  cudaFree(arr);
}

void free_cuda_1i(int *arr)
{
  cudaFree(arr);
}

