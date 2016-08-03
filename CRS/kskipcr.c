#include "kskipcr.h"

void KSKIPCR_init(double *Ar, double *Ap, double *delta, double *eta, double *zeta, double *rvec, double *pvec, double *Av, double *xvec, const int N, const int kskip)
{
  int size = N*(2*kskip+1);
  vec_init(Ar, 0.0, size);
  vec_init(Ap, 0.0, size);
  vec_init(delta, 0.0, size);
  vec_init(eta, 0.0, size);
  vec_init(zeta, 0.0, size);
  vec_init(rvec, 0.0, N);
  vec_init(pvec, 0.0, N);
  vec_init(Av, 0.0, N);
  vec_init(xvec, 0.0, N);
}

int KSKIPCR_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, const struct Parameter *para, const int N, const int NNZ, const bool f_isinner)
{
  int nloop, iloop, jloop;

  double *Ar, *Ap, *delta, *eta, *zeta;
  double *rvec, *pvec, *Av, *x_0;
  double rnorm, bnorm, alpha, beta;
  double error = 0.0;

  bool exit_flag = 2;

  double t_error;

  FILE *p_x=NULL, *p_his=NULL;

  int i_max = 0;
  int i_kskip = 0;
  double d_eps = 0.0;

  bool f_cuda = para->f_cuda;
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
    d_eps = para->d_outer_eps;
  }else{
    i_max = para->i_inner_maxloop;
    i_kskip = para->i_inner_kskip;
    d_eps = para->d_inner_eps;
  }

  if(f_cuda)
  {
    error_log("not done yet");
  }else{
    Ar=malloc_1d(N*(2*i_kskip+1));
    Ap=malloc_1d(N*(2*i_kskip+1));
    delta=malloc_1d(N*(2*i_kskip+1));
    eta=malloc_1d(N*(2*i_kskip+1));
    zeta=malloc_1d(N*(2*i_kskip+1));

    rvec=malloc_1d(N);
    pvec=malloc_1d(N);
    Av=malloc_1d(N);
    x_0=malloc_1d(N);
  }

  //
  //start here
  //

  //init
  KSKIPCR_init(Ar, Ap, delta, eta, zeta, rvec, pvec, Av, xvec, N, i_kskip);

  //backup xvec
  vec_copy(x_0, xvec, N);

  //Ax
  if(f_cuda)
  {
    error_log("not done yet");
  }else{
    MV_mult_CSR(Av, val, col, ptr, xvec, N);
  }

  //r=b-Ax
  vec_sub(rvec, bvec, Av, N);

  //p=r
  vec_copy(pvec, rvec, N);

  //2 norm
  bnorm=norm_2_d(bvec, N);

  for(nloop=0; nloop<i_max; nloop+=(i_kskip+1))
  {
    //2 norm for r
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

    //Ar-> Ar^2k+1
    //Ap-> Ap^2k+1
    cal_arap_kskipcr_d(Ar, Ap, val, col, ptr, rvec, pvec, N, i_kskip);

    //delta=(r,Ar)
    //eta=(A1p,Ap)
    //zeta=(r,Ap)
    cal_deltaetazeta_kskipcr_d(delta, eta, zeta, Ar, Ap, rvec, N, i_kskip);

    for(iloop=nloop; iloop<=nloop+i_kskip; iloop++)
    {
      //alpha = delta_1/eta_1
      alpha=delta[0]/eta[0];

      //beta (delta_1 - 2*alpha*zeta_2 + alpha^2*eta[2])/delta_1
      beta=(delta[0] - 2*alpha*zeta[1] + alpha*alpha*eta[1]) / delta[0];

      //update delta eta zeta
      for(jloop=0; jloop<2*i_kskip-2*(iloop-nloop); jloop++)
      {
        double delta2=0.0;
        delta[jloop] = delta[jloop] - 2*alpha*zeta[jloop+1] + alpha*alpha*eta[jloop+1];
        delta2=delta[jloop+1] - 2*alpha*zeta[jloop+2] + alpha*alpha*eta[jloop+2];
        eta[jloop] = delta2 + 2*beta*(zeta[jloop+1]-alpha*eta[jloop+1]) + beta*beta*eta[jloop];
        zeta[jloop] = delta[jloop] - alpha*zeta[jloop+1] - alpha*(zeta[jloop+1]-alpha*eta[jloop+1]) + beta*zeta[jloop] - alpha*beta*eta[jloop];
      }

      //new Ap
      MV_mult_CSR(Av, val, col, ptr, pvec, N);

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
    free_1d(Ar);
    free_1d(Ap);
    free_1d(delta);
    free_1d(eta);
    free_1d(zeta);
    free_1d(rvec);
    free_1d(pvec);
    free_1d(Av);
    free_1d(x_0);
  }
  if(f_isinner)
  {
    fclose(p_x);
    fclose(p_his);
  }
  return exit_flag;
}
