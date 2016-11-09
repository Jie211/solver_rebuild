#include "cg.h"

void CG_Init(double *rvec, double *pvec, double *Av, double *xvec, const int N)
{
  vec_init(rvec, 0.0, N);
  vec_init(pvec, 0.0, N);
  vec_init(Av, 0.0, N);
  vec_init(xvec, 0.0, N);
}

int CG_CRS(double *val, int *col, int *ptr, double *Tval, int *Tcol, int *Tptr, double *bvec, double *xvec, const struct Parameter *para, const int N, const int NNZ, const bool f_isinner)
{
  int loop;

  double *rvec, *pvec, *Av, *x_0, dot, error=0.0;
  double alpha, beta, bnorm, rnorm;
  double rr, rr2;

  int i_max=0;
  double d_eps=0.0;
  bool f_isVP=para->isVP;
  bool f_verbose=para->f_verbose;
  bool f_cuda=para->f_cuda;

  int exit_flag=2;

  double t_error = 0.0;

  FILE *p_x=NULL, *p_his=NULL;

  if(!f_isinner)
  {
    p_x=file_init("./output/CG_x.txt", "w");
    p_his=file_init("./output/CG_his.txt", "w");
  }

  if(f_cuda)
  {
    error_log("not done yet");
  }else{
    Av = malloc_1d(N);
    rvec = malloc_1d(N);
    pvec = malloc_1d(N);
    x_0 = malloc_1d(N);
  }

  //
  // start here
  //

  //init
  CG_Init(rvec, pvec, Av, xvec, N);

  if(f_isVP && f_isinner)
  {
    i_max = para->i_inner_maxloop;
    d_eps = para->d_inner_eps;
  }else{
    i_max = para->i_outer_maxloop;
    d_eps = para->d_outer_eps;
  }

  //backup xvec
  vec_copy(x_0, xvec, N);

  //b 2norm
  bnorm = norm_2_d(bvec, N);

  //Ax
  if(f_cuda)
  {
    error_log("not done yet");
  }else{
    MV_mult_CSR(Av, val, col, ptr, xvec, N);
  }
  
  //r = b - Ax
  vec_sub(rvec, bvec, Av, N);

  //p = r
  vec_copy(pvec, rvec, N);

  //r dot
  if(f_cuda)
  {
    error_log("not done yet");
  }else{
    rr = dot_d(rvec, rvec, N);
  }
  
  for(loop=0;loop<i_max;loop++)
  {
    rnorm = norm_2_d(rvec, N);
    error = rnorm / bnorm;
    if(!f_isinner)
    {
      if(f_verbose)
      {
        printf("%d %.12e\n", loop+1, error);
      }
      fprintf(p_his, "%d %.12e\n", loop+1, error);
    }
    if(error <= d_eps)
    {
/* #ifdef EBUG */
/*       normal_log("CG convergence"); */
/* #endif */
      exit_flag = 1;
      break;
    }

    //Ap
    if(f_cuda)
    {
      error_log("not done yet");
    }else{
      MV_mult_CSR(Av, val, col, ptr, pvec, N);
    }

    //alpha = (r,r)/(p,ap)
    if(f_cuda)
    {
      error_log("not done yet");
    }else{
      dot = dot_d(pvec, Av, N);
    }

    alpha = rr / dot;

    //x = x + alpha * pvec
    scalar_xpy_d(xvec, alpha, pvec, xvec, N);

    //r = r - alpha * Ap
    scalar_xpy_d(rvec, -alpha, Av, rvec, N);

    //rr2 dot
    if(f_cuda)
    {
      error_log("not done yet");
    }else{
      rr2 = dot_d(rvec, rvec, N);
    }

    beta = rr2 / rr;

    rr = rr2;

    //p = r + beta * p
    scalar_xpy_d(pvec, beta, pvec, rvec, N);
  }

  if(!f_isinner)
  {
    file_print(p_x, xvec, N);
    fclose(p_x);
    fclose(p_his);
    t_error = error_d_CRS(val, col, ptr, bvec, xvec, x_0, N);
    printf("|b-ax|2/|b|2=%.1f\n", t_error);
    printf("loop=%d\n", loop+1);
  }else{
    /* fclose(p_x); */
    /* fclose(p_his); */
  }
  if(f_isinner && f_verbose)
  {
    printf("Inner %d %.12e\n", loop+1, error);
  }

  if(f_cuda)
  {
    error_log("not done yet");
  }else{
    free_1d(rvec);
    free_1d(pvec);
    free_1d(Av);
    free_1d(x_0);
  }
  
  return exit_flag;
}
