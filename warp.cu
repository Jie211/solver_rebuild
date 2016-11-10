#include "warp.h"

__global__ void cuda_kernel()
{

}

void cuda_kernel_warp(){
  cuda_kernel<<<1, 1>>>();
}
