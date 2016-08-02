#include "cr.h"
void CR_init(double *rvec, double *pvec, double *qvec, double *svec, double *xvec, const int N)
{
  vec_init(rvec, 0.0, N);
  vec_init(pvec, 0.0, N);
  vec_init(qvec, 0.0, N);
  vec_init(svec, 0.0, N);
  vec_init(xvec, 0.0, N);
}

int CR_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, const struct Parameter *para, const int N, const int NNZ, const bool isinner)
{
  int loop;

  double *rvec, *pvec, *qvec, *svec, *x_0, error = 0.0;
  double alpha, beta, bnorm, rnorm;
  double rs, rs2;

  int  exit_flag = 2;
  int i_max = 0;
  double d_eps = 0.0;

  double t_error = 0.0;

  FILE *p_x=NULL, *p_his=NULL;

  if(!isinner)
  {
    p_x=file_init("./output/CR_x.txt", "w");
    p_his=file_init("./output/CR_his.txt", "w");
  }

  if(para->f_cuda != true)
  {
    rvec = malloc_1d(N);
    pvec = malloc_1d(N);
    qvec = malloc_1d(N);
    svec = malloc_1d(N);
    x_0 = malloc_1d(N);
  }else{
    error_log("not done yet");
  }

  //
  //start here
  //

  //int
  CR_init(rvec, pvec, qvec, svec, xvec, N);

  if(para->isVP == true && isinner == true)
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
  if(para->f_cuda!=true)
  {
    MV_mult_CSR(qvec, val, col, ptr, xvec, N);
  }else{
    error_log("not done yet");
  }

  //r=b-Ax
  vec_sub(rvec, bvec, qvec, N);

  //p=r
  vec_copy(pvec, rvec, N);

  //q=Ap
  if(para->f_cuda!=true)
  {
    MV_mult_CSR(qvec, val, col, ptr, pvec, N);
  }

  //s=q
  vec_copy(svec, qvec, N);

  //(r, s)
  if(para->f_cuda!=true)
  {
    rs = dot_d(rvec, svec, N);
  }

  for(loop=0;loop<i_max;loop++)
  {
    rnorm = norm_2_d(rvec, N);
    error = rnorm / bnorm;
    if(!isinner)
    {
      if(para->f_verbose)
      {
        printf("%d %.12e\n", loop+1, error);
      }
      fprintf(p_his, "%d %.12e\n", loop+1, error);
    }
    if(error <= d_eps)
    {
#ifdef EBUG
      normal_log("CR convergence");
#endif
      exit_flag = 1;
      break;
    }

    //alpha = (r,s)/(q,q)
    alpha = rs / dot_d(qvec, qvec, N);

    //x=x+alpha*pvec
    scalar_xpy_d(xvec, alpha, pvec, xvec, N);

    //r=r-alpha*qvec
    scalar_xpy_d(rvec, -alpha, qvec, rvec, N);

    //s=Ar
    if(para->f_cuda!=true)
    {
      MV_mult_CSR(svec, val, col, ptr, rvec, N);
    }

    //(r, s)
    rs2 = dot_d(rvec, svec, N);

    //beta=(r_new,s_new)/(r,s)
    beta = rs2/rs;

    rs = rs2;

    //p=r+beta*p
    scalar_xpy_d(pvec, beta, pvec, rvec, N);

    //q=s+beta*q
    scalar_xpy_d(qvec, beta, qvec, svec, N);
  }

  if(!isinner)
  {
    file_print(p_x, xvec, N);
    fclose(p_x);
    t_error = error_d_CRS(val, col, ptr, bvec, xvec, x_0, N);
    printf("|b-ax|2/|b|2=%.1f\n", t_error);
    printf("loop=%d\n", loop+1);
  }
  if(isinner && para->f_verbose)
  {
    printf("Inner %d %.12e\n", loop+1, error);
  }
  if(para->f_cuda != true)
  {
    free_1d(rvec);
    free_1d(pvec);
    free_1d(qvec);
    free_1d(x_0);
  }
  if(isinner)
  {
    fclose(p_x);
    fclose(p_his);
  }
  return exit_flag;
}