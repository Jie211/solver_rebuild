#include "kskipcg.h"

void KSKIPCG_init(double **Ap, double **Ar, double *delta, double *eta, double *zeta, double *rvec, double *pvec, double *Av, double *xvec, const int N, const int kskip)
{
  vec_init_2(Ar, 0.0, N, (2*kskip+1));
  vec_init_2(Ap, 0.0, N, (2*kskip+2));
  vec_init(delta, 0.0, 2*kskip);
  vec_init(eta, 0.0, 2*kskip+1);
  vec_init(zeta, 0.0, 2*kskip+2);
  vec_init(rvec, 0.0, N);
  vec_init(pvec, 0.0, N);
  vec_init(Av, 0.0, N);
  vec_init(xvec, 0.0, N);
}

int KSKIPCG_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, const struct Parameter *para, const int N, const int NNZ, const bool f_isinner)
{
  int nloop, iloop, jloop;
  /* int i, ii; */

  double *delta, *eta, *zeta;
  double *rvec, *pvec, *Av, *x_0;
  double **Ap, **Ar;

  double rnorm, bnorm, alpha, beta, gamma;

  int exit_flag = 2;
  double error = 0.0;

  double t_error;

  FILE *p_x=NULL, *p_his=NULL;

  /* double dot=0.0, tmp1=0.0, tmp2=0.0, tmp3=0.0; */

  int i_max = 0;
  double d_eps = 0.0;
  int i_kskip = 0;
  bool f_cuda = para->f_cuda;
  int i_fix;
  bool f_verbose = para->f_verbose;

  if(!f_isinner)
  {
    p_x=file_init("./output/KSKIPCG_x.txt", "w");
    p_his=file_init("./output/KSKIPCG_his.txt", "w");
  }

  if(!f_isinner)
  {
    i_max = para->i_outer_maxloop;
    i_kskip = para->i_outer_kskip;
    i_fix = para->i_outer_fix;
    d_eps = para->d_outer_eps;
  }else{
    i_max = para->i_inner_maxloop;
    i_kskip = para->i_inner_kskip;
    i_fix = para->i_inner_fix;
    d_eps = para->d_inner_eps;
  }

  if(f_cuda)
  {
    error_log("not done yet");
  }else{
    Ar=malloc_2d(N, 2*i_kskip+1);
    Ap=malloc_2d(N, 2*i_kskip+2);
    
    delta=malloc_1d(2*i_kskip);
    eta=malloc_1d(2*i_kskip+1);
    zeta=malloc_1d(2*i_kskip+2);

    rvec=malloc_1d(N);
    pvec=malloc_1d(N);
    Av=malloc_1d(N);
    x_0=malloc_1d(N);
  }
  
  //
  //start here
  //

  //init
  KSKIPCG_init(Ap, Ar, delta, eta, zeta, rvec, pvec, Av, xvec, N, i_kskip);

  //backup xvec
  vec_copy(x_0, xvec, N);

  //Av
  if(f_cuda)
  {
    error_log("not done yet 1");
  }else{
    MV_mult_CSR(Av, val, col, ptr, xvec, N);
  }

  //r=b-Ax
  vec_sub(rvec, bvec, Av, N);

  //p=r
  vec_copy(pvec, rvec, N);

  //b 2norm
  bnorm = norm_2_d(bvec, N);

  for(nloop=0;nloop<i_max;nloop+=(i_kskip+1))
  {
    rnorm=norm_2_d(rvec, N);
    error=rnorm/bnorm;
    if(!f_isinner)
    {
      if(f_verbose)
      {
        printf("%d %.12e\n", nloop+1, error);
      }
      fprintf(p_his, "%d %.12e\n", nloop+1, error);
    }
    if(error<=d_eps)
    {
      exit_flag = 1;
      break;
    }

    //Ar-> Ar^2k
    //Ap-> Ap^2k+2
    if(f_cuda)
    {
      error_log("not done yet 2");
    }else{
      cal_arap_kskipcg_d(Ar, Ap, val, col, ptr, rvec, pvec, N, i_kskip);
    }

    //gamma=(r,r)
    if(f_cuda)
    {
      error_log("not done yet 3");
    }else{
      gamma = dot_d(rvec, rvec, N);
    }

    //delta=(r,Ar)
    //eta=(r,Ap)
    //zeta=(p,Ap)
    if(f_cuda)
    {
      error_log("not done yet 4");
    }else{
      cal_deltaetazeta_kskipcg_d(delta, eta, zeta, Ar, Ap, rvec, pvec, N, i_kskip);
    }

    for(iloop=nloop;iloop<=nloop+i_kskip;iloop++)
    {
      //alpha = gamma/zeta_1
      alpha=gamma/zeta[0];
      //beta = (alpha*zeta_2/zeta_1) -1
      beta=alpha*zeta[1]/zeta[0]-1.0;
      if(i_fix==1)
      {
        gamma=beta*gamma;
      }else if(i_fix==2)
      {
        double tmp0=gamma-alpha*eta[0];
        double tmp1=eta[0]-alpha*zeta[1];
        gamma=tmp0-alpha*tmp1;
      }

      //update delta eta zeta
      for(jloop=0; jloop<2*i_kskip-2*(iloop-nloop);jloop++){
        delta[jloop] = delta[jloop] - 2*alpha*eta[jloop+1] + alpha*alpha*eta[jloop+2];
        double eta_old=eta[jloop];
        eta[jloop] = delta[jloop] + beta*zeta[jloop+1] - alpha*beta*zeta[jloop+1];
        zeta[jloop] = eta[jloop+1] + beta*eta_old + beta*beta*zeta[jloop] - alpha*beta*zeta[jloop+1];
      }

      if(f_cuda)
      {
        error_log("not done yet 5");
      }else{
        MV_mult_CSR(Av, val, col, ptr, pvec, N);
      }

      //x=x+alpha*p
      scalar_xpy_d(xvec, alpha, pvec, xvec, N);

      //r=r-alpha*Ap
      scalar_xpy_d(rvec, -alpha, Av, rvec, N);

      //p=r+beta*p
      scalar_xpy_d(pvec, beta, pvec, rvec, N);
    }
  }
  if(!f_isinner)
  {
    file_print(p_x, xvec, N);
    fclose(p_x);
    t_error = error_d_CRS(val, col, ptr, bvec, xvec, x_0, N);
    printf("|b-ax|2/|b|2=%.1f\n", t_error);
    printf("loop=%d\n", nloop+1);
  }
  if(f_isinner && f_verbose)
  {
    printf("Inner %d %.12e\n", nloop+1, error);
  }
  if(f_cuda)
  {
    error_log("not done yet 6");
  }else{
    free_1d(delta);
    free_1d(eta);
    free_1d(zeta);
    free_1d(rvec);
    free_1d(pvec);
    free_1d(Av);
    free_1d(x_0);
    free_2d(Ap, 2*i_kskip+2);
    free_2d(Ar, 2*i_kskip+1);
  }
  if(f_isinner)
  {
    fclose(p_x);
    fclose(p_his);
  }
  return exit_flag;
}
