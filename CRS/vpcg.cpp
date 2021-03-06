#include "vpcg.h"

void VPCG_Init(double *rvec, double *pvec, double *Av, double *zvec, double *xvec, const int N)
{
  vec_init(rvec, 0.0, N);
  vec_init(pvec, 0.0, N);
  vec_init(Av, 0.0, N);
  vec_init(zvec, 0.0, N);
  vec_init(xvec, 0.0, N);
}

int VPCG_CRS(double *val, int *col, int *ptr, double *Tval, int *Tcol, int *Tptr, double *bvec, double *xvec, struct Parameter *para, const int N, const int NNZ, const bool f_isinner)
{
  int loop;

  double *rvec, *pvec, *Av, *x_0, dot, error=0.0;
  double *zvec;
  double alpha, beta, bnorm, rnorm;
  double rz, rz2;

  int i_max=0;
  double d_eps=0.0;
  /* bool f_isVP=para->isVP; */
  bool f_verbose=para->f_verbose;
  bool f_cuda=para->f_cuda;

  int exit_flag=2;

  double t_error = 0.0;

  FILE *p_x=NULL, *p_his=NULL;

  int get_error;

  /* if(!f_isinner) */
  /* { */
  p_x=file_init("./output/VPCG_x.txt", "w");
  p_his=file_init("./output/VPCG_his.txt", "w");
  /* } */

  if(f_cuda)
  {
    error_log((char*)"not done yet");
  }else{
    Av = malloc_1d(N);
    rvec = malloc_1d(N);
    pvec = malloc_1d(N);
    zvec = malloc_1d(N);
    x_0 = malloc_1d(N);
  }

  //
  // start here
  //

  //init
  VPCG_Init(rvec, pvec, Av, zvec, xvec, N);

  /* if(f_isVP && f_isinner) */
  /* { */
  /*   i_max = para->i_inner_maxloop; */
  /*   d_eps = para->d_inner_eps; */
  /* }else{ */
  /*   i_max = para->i_outer_maxloop; */
  /*   d_eps = para->d_outer_eps; */
  /* } */

  i_max = para->i_outer_maxloop;
  d_eps = para->d_outer_eps;


  //backup xvec
  vec_copy(x_0, xvec, N);

  //b 2norm
  bnorm = norm_2_d(bvec, N);

  //Ax
  if(f_cuda)
  {
    error_log((char*)"Cuda not done yet");
  }else{
    MV_mult_CSR(Av, val, col, ptr, xvec, N);
  }
  
  //r = b - Ax
  vec_sub(rvec, bvec, Av, N);

  //---------->

  // Az=r
  get_error = inner_selecter(para, rvec, zvec, val, col, ptr, Tval, Tcol, Tptr, N, NNZ);
  if(get_error==-1)
  {
    error_log((char*)"errror in vpcg - inner_selecter");
  }

  //p = z
  vec_copy(pvec, zvec, N);

  //r, z
  if(f_cuda)
  {
    error_log((char*)"not done yet");
  }else{
    rz = dot_d(rvec, zvec, N);
  }
  
  for(loop=0;loop<i_max;loop++)
  {
    rnorm = norm_2_d(rvec, N);
    error = rnorm / bnorm;
    if(f_verbose)
    {
      printf("Outer %d %.12e\n", loop+1, error);
    }
    fprintf(p_his, "%d %.12e\n", loop+1, error);
    if(error <= d_eps)
    {
      exit_flag = 1;
      break;
    }

    //Ap
    if(f_cuda)
    {
      error_log((char*)"not done yet");
    }else{
      MV_mult_CSR(Av, val, col, ptr, pvec, N);
    }

    //alpha = (r,r)/(p,ap)
    if(f_cuda)
    {
      error_log((char*)"not done yet");
    }else{
      dot = dot_d(pvec, Av, N);
    }

    alpha = rz / dot;

    //x = x + alpha * pvec
    scalar_xpy_d(xvec, alpha, pvec, xvec, N);

    //r = r - alpha * Ap
    scalar_xpy_d(rvec, -alpha, Av, rvec, N);

    vec_init(zvec, 0.0, N);

    get_error = inner_selecter(para, rvec, zvec, val, col, ptr, Tval, Tcol, Tptr, N, NNZ);
    if(get_error==-1)
    {
      error_log((char*)"errror in vpcg - inner_selecter");
    }


    //zr2 dot
    if(f_cuda)
    {
      error_log((char*)"not done yet");
    }else{
      rz2 = dot_d(rvec, zvec, N);
    }

    beta = rz2 / rz;

    rz = rz2;

    //p = z + beta * p
    scalar_xpy_d(pvec, beta, pvec, zvec, N);
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
    free_1d(Av);
    free_1d(x_0);
    free_1d(zvec);
  }
  
  return exit_flag;
}
