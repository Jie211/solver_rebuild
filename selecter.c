#include "selecter.h"

int outer_selecter(struct Parameter *para, double *bvec, double *xvec, double *val, int *col, int *ptr, const int N, const int NNZ)
{
  int handle = 0;
  bool isinner = false;

  if(para->isVP!=true)
  {
    if(para->c_outer_solver == CG)
    {
      handle = CG_CRS(val, col, ptr, bvec, xvec, para, N, NNZ, isinner);
    }else if(para->c_outer_solver == CR)
    {
      handle = CR_CRS(val, col, ptr, bvec, xvec, para, N, NNZ, isinner);
    }else if(para->c_outer_solver == GCR)
    {
      handle = GCR_CRS(val, col, ptr, bvec, xvec, para, N, NNZ, isinner);
    }else if(para->c_outer_solver == GMRES)
    {
      handle = GMRES_CRS(val, col, ptr, bvec, xvec, para, N, NNZ, isinner);
    }
    else{
      error_log("not define now");
      return -1;
    }
  }else{
    error_log("not define now");
    return -1;
  }
  if(handle==1)
  {
    normal_log("OuterSolver convergence");
  }else{
    normal_log("OuterSolver NOT convergence");
  }
  return 0;
}
