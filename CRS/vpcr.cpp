#include "vpcr.h"

void VPCR_Init(double *rvec, double *pvec, double *zvec, double *Av, double *Ap, double *xvec, const int N)
{
  vec_init(rvec, 0.0, N);
  vec_init(pvec, 0.0, N);
  vec_init(zvec, 0.0, N);
  vec_init(Av, 0.0, N);
  vec_init(Ap, 0.0, N);
  vec_init(xvec, 0.0, N);
}

int VPCR_CRS(double *val, int *col, int *ptr, double *Tval, int *Tcol, int *Tptr, double *bvec, double *xvec, struct Parameter *para, const int N, const int NNZ, const bool f_isinner)
{
  int loop;

  double *rvec, *pvec, *zvec, *Av, *Ap, *x_0, error = 0.0;
  double alpha, beta, bnorm, rnorm;
  double zaz, zaz2;
  double tmp;

  int i_max = 0;
  double d_eps = 0.0;
  bool f_verbose = para->f_verbose;
  bool f_cuda = para->f_cuda;
  
  int exit_flag = 2;

  double t_error = 0.0;
  FILE *p_x = NULL, *p_his = NULL;

  int get_error;

  p_x = file_init("./output/VPCR_x.txt", "w");
  p_his = file_init("./output/VPCR_his.txt", "w");

  if(f_cuda)
  {
    error_log((char*)"cuda not done yet");
  }else{
    rvec = malloc_1d(N);
    pvec = malloc_1d(N);
    zvec = malloc_1d(N);
    Av = malloc_1d(N);
    Ap = malloc_1d(N);
    x_0 = malloc_1d(N);
  }

  //
  //start here
  //

  //init
  VPCR_Init(rvec, pvec, zvec, Av, Ap, xvec, N);

  i_max = para->i_outer_maxloop;
  d_eps = para->d_outer_eps;

  //backup xvec
  vec_copy(x_0, xvec, N);

  //b 2norm
  bnorm = norm_2_d(bvec, N);

  if(f_cuda){
    error_log((char*)"Cuda not done yet");
  }else{
    MV_mult_CSR(Av, val, col, ptr, xvec, N);
  }

  //r=b-Ax
  vec_sub(rvec, bvec, Av, N);

  // Az=r
  get_error = inner_selecter(para, rvec, zvec, val, col, ptr, Tval, Tcol, Tptr, N, NNZ);
  if(get_error==-1)
  {
    error_log((char*)"error in vpcr - inner_selecter");
  }

  //p=z
  vec_copy(pvec, zvec, N);

  //Az(Av)
  if(f_cuda){
    error_log((char*)"Cuda not done yet");
  }else{
    MV_mult_CSR(Av, val, col, ptr, zvec, N);
  }


  //Ap=Az
  vec_copy(Ap, Av, N);

  zaz=dot_d(zvec, Av, N);

  for(loop=0; loop<i_max; loop++)
  {
    rnorm = norm_1_d(rvec, N);
    error = rnorm / bnorm;
    if(f_verbose)
    {
      printf("Outer %d %.12e\n", loop+1, error);
    }
    fprintf(p_his, "%d %.12e\n", loop+1, error);;
    if(error <= d_eps)
    {
      exit_flag = 1;
      break;
    }

    //alpha = (z,Az)/(Ap,Ap)
    if(f_cuda)
    {
      error_log((char*)"cuda not done yet");
    }else{
      tmp = dot_d(Ap, Ap, N);
    }
    alpha = zaz / tmp;

    //x=alpha*pvec+x
    scalar_xpy_d(xvec, alpha, pvec, xvec, N);

    //r=-alpha*Ap+r
    scalar_xpy_d(rvec, -alpha, Ap, rvec, N);

    //init zvec
    vec_init(zvec, 0.0, N);

    //Az = r
    get_error = inner_selecter(para, rvec, zvec, val, col, ptr, Tval, Tcol, Tptr, N, NNZ);
    if(get_error==-1)
    {
      error_log((char*)"error in vpcg - inner_selecter");
    }

    //Az
    if(f_cuda){
      error_log((char*)"Cuda not done yet");
    }else{
      MV_mult_CSR(Av, val, col, ptr, zvec, N);
    }

    //(z,Az)
    if(f_cuda){
      error_log((char*)"cuda done yet");
    }else{
      zaz2=dot_d(zvec, Av, N);
    }

    beta = zaz2/zaz;

    zaz = zaz2;

    //p=beta*p+z
    scalar_xpy_d(pvec, beta, pvec, zvec, N);

    //Ap=beta*Ap+Az
    scalar_xpy_d(Ap, beta, Ap, Av, N);
  }

  file_print(p_x, xvec, N);
  fclose(p_x);
  fclose(p_his);
  t_error = error_d_CRS(val, col, ptr, bvec, xvec, x_0, N);
  printf("|b-ax|2/|b|2=%.1f\n", t_error);
  printf("loop=%d\n", loop+1);
  if(f_cuda)
  {
    error_log((char*)"not done yet");
  }else{
    free_1d(rvec);
    free_1d(pvec);
    free_1d(zvec);
    free_1d(Av);
    free_1d(Ap);
    free_1d(x_0);
  }
  
  return exit_flag;
}

