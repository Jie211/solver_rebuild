#include "selecter.h"

int outer_selecter(struct Parameter *para, double *bvec, double *xvec, double *val, int *col, int *ptr, double *Tval, int *Tcol, int *Tptr, const int N, const int NNZ)
{
  int handle_out = 0;
  bool isinner = false;

  if(para->isVP==false)
  {
    if(para->c_outer_solver == CG)
    {
      handle_out = CG_CRS(val, col, ptr, Tval, Tcol, Tptr, bvec, xvec, para, N, NNZ, isinner);
    }else if(para->c_outer_solver == CR)
    {
      handle_out = CR_CRS(val, col, ptr, Tval, Tcol, Tptr, bvec, xvec, para, N, NNZ, isinner);
    }else if(para->c_outer_solver == GCR)
    {
      handle_out = GCR_CRS(val, col, ptr, Tval, Tcol, Tptr, bvec, xvec, para, N, NNZ, isinner);
    }else if(para->c_outer_solver == GMRES)
    {
      handle_out = GMRES_CRS(val, col, ptr, Tval, Tcol, Tptr, bvec, xvec, para, N, NNZ, isinner);
    }else if(para->c_outer_solver == KSKIPCG)
    {
      handle_out = KSKIPCG_CRS(val, col, ptr, Tval, Tcol, Tptr, bvec, xvec, para, N, NNZ, isinner);
    }else if(para->c_outer_solver == KSKIPCR)
    {
      warning_log((char*)"kskipcr have some bug");
    }else if(para->c_outer_solver == BICG)
    {
      handle_out = BICG_CRS(val, col, ptr, Tval, Tcol, Tptr, bvec, xvec, para, N, NNZ, isinner);
    }else if(para->c_outer_solver == KSKIPBICG)
    {
      handle_out = KSKIPBICG_CRS(val, col, ptr, Tval, Tcol, Tptr, bvec, xvec, para, N, NNZ, isinner);
    }
    else{
      error_log((char*)"not define now");
      return -1;
    }
  }else if(para->isVP==true)
  {
    if(para->c_outer_solver == VPCG)
    {
      handle_out = VPCG_CRS(val, col, ptr, Tval, Tcol, Tptr, bvec, xvec, para, N, NNZ, isinner);
    }else if(para->c_outer_solver == VPCR)
    {
      handle_out = VPCR_CRS(val, col, ptr, Tval, Tcol, Tptr, bvec, xvec, para, N, NNZ, isinner);
    }else if(para->c_outer_solver == VPGCR)
    {
      handle_out = VPGCR_CRS(val, col, ptr, Tval, Tcol, Tptr, bvec, xvec, para, N, NNZ, isinner);
    }else if(para->c_outer_solver == VPGMRES)
    {
      handle_out = VPGMRES_CRS(val, col, ptr, Tval, Tcol, Tptr, bvec, xvec, para, N, NNZ, isinner);
    }
    else{
      error_log((char*)"not define now");
      return -1;
    }
  }
  else{
    error_log((char*)"not define now");
    return -1;
  }
  if(handle_out==1)
  {
    normal_log((char*)"OuterSolver convergence");
  }else if(handle_out==-1)
  {
    error_log((char*)"error in outer_selecter");
  }else{
    normal_log((char*)"OuterSolver NOT convergence");
  }
  return 0;
}

int inner_selecter(struct Parameter *para, double *bvec, double *xvec, double *val, int *col, int *ptr, double *Tval, int *Tcol, int *Tptr, int N, int NNZ)
{
  int handle_in = 0;
  bool isinner = true;

  if(para->isVP==true)
  {
    if(para->c_inner_solver == CG)
    {
      handle_in = CG_CRS(val, col, ptr, Tval, Tcol, Tptr, bvec, xvec, para, N, NNZ, isinner);
    }else if(para->c_inner_solver == CR)
    {
      handle_in = CR_CRS(val, col, ptr, Tval, Tcol, Tptr, bvec, xvec, para, N, NNZ, isinner);
    }else if(para->c_inner_solver == GCR)
    {
      handle_in = GCR_CRS(val, col, ptr, Tval, Tcol, Tptr, bvec, xvec, para, N, NNZ, isinner);
    }else if(para->c_inner_solver == GMRES)
    {
      handle_in = GMRES_CRS(val, col, ptr, Tval, Tcol, Tptr, bvec, xvec, para, N, NNZ, isinner);
    }else if(para->c_inner_solver == KSKIPCG)
    {
      handle_in = KSKIPCG_CRS(val, col, ptr, Tval, Tcol, Tptr, bvec, xvec, para, N, NNZ, isinner);
    }else if(para->c_inner_solver == KSKIPCR)
    {
      warning_log((char*)"kskipcr have some bug");
    }else if(para->c_inner_solver == BICG){
      handle_in = BICG_CRS(val, col, ptr, Tval, Tcol, Tptr, bvec, xvec, para, N, NNZ, isinner);
    }else if(para->c_inner_solver == KSKIPBICG){
      handle_in = KSKIPBICG_CRS(val, col, ptr, Tval, Tcol, Tptr, bvec, xvec, para, N, NNZ, isinner);
    }
    else{
      error_log((char*)"not define now");
      return -1;
    }
  }else{
    error_log((char*)"error");
    return -1;
  }
  if(handle_in==1)
  {
    /* normal_log((char*)"InnerSolver convergence"); */
  }else if(handle_in==2){
    /* normal_log((char*)"InnerSolver NOT convergence"); */
  }else{
    error_log((char*)"error in inner_selecter");
    return -1;
  }
  return 0;
}
