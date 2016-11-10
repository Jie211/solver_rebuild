#include "start.h"

int csr_start(int argc, char *argv[])
{

  int error;
  struct Parameter para;
  
  int N, NNZ;
  double *bvec, *xvec, *val;
  int *col, *ptr;
  double *Tval;
  int *Tcol, *Tptr;


  init_ver(&para);

  error = get_opt(argc, argv, &para);
  if(error_handle(error, "error in get_cmd")!=0)
    return -1;

  error = check_opt(&para);
  if(error_handle(error, "error in check_cmd")!=0)
    return -1;

  error = find_mat(&para);
  if(error_handle(error, "error in find_mat")!=0)
    return -1;

  show_opt(&para);

  error = set_openmp_thread(para.i_thread);
  if(error_handle(error, "error in set_openmp_thread")!=0)
    return -1;

  error = get_mat_head(&para, &N, &NNZ);
  if(error_handle(error, "error in get_mat_head")!=0)
    return -1;

  if( para.f_cuda == false )
  {
    bvec = malloc_1d(N);
    xvec = malloc_1d(N);

    val = malloc_1d(NNZ);
    col = malloc_1i(NNZ);
    ptr = malloc_1i(N+1);

    Tval=malloc_1d(NNZ);
    Tcol=malloc_1i(NNZ);
    Tptr=malloc_1i(N+1);
  }else{
    error_log("Cuda not done now");
    return -1;
  }


  error = get_mat_data(&para, col, ptr, val, bvec, xvec, N, NNZ);
  if(error_handle(error, "error in get_mat_data")!=0)
    return -1;

  //A^T
  Transpose_d(val, col, ptr, Tval, Tcol, Tptr, N, NNZ);

  error = outer_selecter(&para, bvec, xvec, val, col, ptr, Tval, Tcol, Tptr, N, NNZ);
  if(error_handle(error, "error in outer_selecter")!=0)
    return -1;

  show_opt(&para);

  if( para.f_cuda == false )
  {
    free_1d(bvec);
    free_1d(xvec);

    free_1d(val);
    free_1i(col);
    free_1i(ptr);

    free_1d(Tval);
    free_1i(Tcol);
    free_1i(Tptr);
  }else{
    error_log("Cuda not done now");
    return -1;
  }
  return 0;
}
