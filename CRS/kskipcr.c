#include "kskipcr.h"
void KSKIPCR_Init(double *v1, double *v2, double *v3, double *v4, double *v5, double *v6, double *v7, double *v8, double *v9, int N, int kskip)
{
  vec_init(v1, 0.0, (2*kskip+1)*N);
  vec_init(v2, 0.0, (2*kskip+1)*N);
  vec_init(v3, 0.0, (2*kskip+1)*N);
  vec_init(v4, 0.0, (2*kskip+1)*N);
  vec_init(v5, 0.0, (2*kskip+1)*N);

  vec_init(v6, 0.0, N);
  vec_init(v7, 0.0, N);
  vec_init(v8, 0.0, N);
  vec_init(v9, 0.0, N);
}

int KSKIPCR_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, const struct Parameter *para, const int N, const int NNZ, const bool f_isinner)
{

  int nloop, iloop, jloop;

  double *Ar, *Ap, *delta, *eta, *zeta;
  double *rvec, *pvec, *Av, *x_0;
  double rnorm, bnorm, alpha, beta, error=0.0;
  bool flag=false;

  double t_error;
  FILE *p_x=NULL, *p_his=NULL;

  int i_max = 0;
  int i_kskip = 0;
  double d_eps = 0.0;

  bool f_cuda = para->f_cuda;
  bool f_verbose = para->f_verbose;
  int i, j, tmp, ii, tmp1, tmp2, tmp3;


  if(!f_isinner){
    p_x=file_init("./output/KSKIPCR_x.txt", "w");
    p_his=file_init("./output/KSKIPCR_his.txt", "w");
  }

  if(!f_isinner)
  {
    i_max=para->i_outer_maxloop;
    i_kskip=para->i_outer_kskip;
    d_eps=para->d_outer_eps;
  }

  Ar = (double *)malloc(sizeof(double) * (2*i_kskip+1)*N );
  Ap = (double *)malloc(sizeof(double) * (2*i_kskip+1)*N );
  delta = (double *)malloc(sizeof(double) * (2*i_kskip+1)*N );
  eta = (double *)malloc(sizeof(double) * (2*i_kskip+1)*N );
  zeta = (double *)malloc(sizeof(double) * (2*i_kskip+1)*N );

  rvec = (double *)malloc(sizeof(double) * N );
  pvec = (double *)malloc(sizeof(double) * N );
  Av = (double *)malloc(sizeof(double) * N );
  x_0 = (double *)malloc(sizeof(double) * N );

  for(i=0;i < (2*i_kskip+1)*N ;i++ )
  {
    Ar[i]=0.0;
    Ap[i]=0.0;
    delta[i]=0.0;
    eta[i]=0.0;
    zeta[i]=0.0;
  }
  for(i=0;i<N;i++)
  {
    rvec[i]=0.0;
    pvec[i]=0.0;
    Av[i]=0.0;
    xvec[i]=0.0;
  } 

  for(i=0;i<N;i++)
  {
    x_0[i]=xvec[i];
  }

  MV_mult_CSR(Av, val, col, ptr, xvec, N);

  /* for(i=0;i<N;i++) */
  /* { */
  /*   tmp = 0.0; */
  /*   for(j=ptr[i];j<ptr[i+1];j++) */
  /*   { */
  /*     tmp+=val[j]*xvec[col[j]]; */
  /*   } */
  /*   Av[i] = tmp; */
  /* } */
  

  for(i=0;i<N;i++)
  {
    rvec[i] = bvec[i]-Av[i];
  }

  for(i=0;i<N;i++)
  {
    pvec[i]=rvec[i];
  }

  bnorm = norm_2_d(bvec, N);

  /* tmp = 0.0; */
  /* for(i=0;i<N;i++) */
  /* { */
  /*   tmp+=bvec[i]*bvec[i]; */
  /* } */
  /* bnorm = sqrt(tmp); */

  for(nloop=0;nloop<i_max;nloop+=(i_kskip+1))
  {
    rnorm=norm_2_d(rvec, N);

    /* tmp = 0.0; */
    /* for(i=0;i<N;i++) */
    /* { */
    /*   tmp+=rvec[i]*rvec[i]; */
    /* } */
    /* rnorm = sqrt(tmp); */


    error = rnorm/bnorm;
    printf("%d %.12e\n",  nloop+1, error);

    if(error <= d_eps)
    {
      break;
    }

    cal_arap_kskipcr_d(Ar, Ap, val, col, ptr, rvec, pvec, N, i_kskip);

    /* for(i=0;i<N;i++) */
    /* { */
    /*   tmp1=0.0; */
    /*   tmp2=0.0; */
    /*   for(j=ptr[i];j<ptr[i+1];j++) */
    /*   { */
    /*     tmp1+=val[j]*rvec[col[j]]; */
    /*     tmp2+=val[j]*pvec[col[j]]; */
    /*   } */
    /*   Ar[i]=tmp1; */
    /*   Ap[i]=tmp2; */
    /* } */
    /* for(ii=1;ii<2*i_kskip+1;ii++) */
    /* { */
    /*   for(i=0;i<N;i++) */
    /*   { */
    /*     tmp1=0.0; */
    /*     tmp2=0.0; */
    /*     for(j=ptr[i];j<ptr[i+1];j++) */
    /*     { */
    /*       tmp1+=val[j]*Ar[(ii-1)*N+col[j]]; */
    /*       tmp2+=val[j]*Ap[(ii-1)*N+col[j]]; */
    /*     } */
    /*     Ar[ii*N+i]=tmp1; */
    /*     Ap[ii*N+i]=tmp2; */
    /*   } */
    /* } */

    cal_deltaetazeta_kskipcr_d(delta, eta, zeta, Ar, Ap, rvec, N, i_kskip);

    /* for(i=0;i<2*i_kskip+1;i++) */
    /* { */
    /*   tmp1=0.0; */
    /*   tmp2=0.0; */
    /*   tmp3=0.0; */
    /*   for(j=0;j<N;j++) */
    /*   { */
    /*     tmp1+=rvec[j]*Ar[i*N+j]; */
    /*     tmp2+=Ap[0*N+j]*Ap[i*N+j]; */
    /*     tmp3+=rvec[j]*Ap[i*N+j]; */
    /*   } */
    /*   delta[i]=tmp1; */
    /*   eta[i]=tmp2; */
    /*   zeta[i]=tmp3; */
    /* } */

    for(iloop=nloop;iloop<=nloop+i_kskip;iloop++)
    {
      alpha = delta[0]/eta[0];

      beta=(delta[0] - 2*alpha*zeta[1]+alpha*alpha*eta[1])/delta[0];

      for(jloop=0;jloop<2*i_kskip-2*(iloop-nloop);jloop++)
      {
        double delta2 = 0.0;
        delta[jloop] = delta[jloop] - 2*alpha*zeta[jloop+1]+alpha*alpha*eta[jloop+1];
        delta2=delta[jloop+1] - 2*alpha*zeta[jloop+2] + alpha*alpha*eta[jloop+2];
        eta[jloop] = delta2 + 2*beta*(zeta[jloop+1]-alpha*eta[jloop+1])+beta*beta*eta[jloop];
        zeta[jloop] = delta[jloop] - alpha*zeta[jloop+1]-alpha*(zeta[jloop+1]-alpha*eta[jloop+1]) + beta*zeta[jloop] - alpha*beta*eta[jloop];
      }

      MV_mult_CSR(Av, val, col, ptr, pvec, N);

      /* for(i=0;i<N;i++) */
      /* { */
      /*   tmp = 0.0; */
      /*   for(j=ptr[i];j<ptr[i+1];j++) */
      /*   { */
      /*     tmp+=val[j]*pvec[col[j]]; */
      /*   } */
      /*   Av[i] = tmp; */
      /* } */

      scalar_xpy_d(xvec, alpha, pvec, xvec, N);

      /* tmp=0.0; */
      /* for(i=0;i<N;i++) */
      /* { */
      /*   tmp=xvec[i]; */
      /*   xvec[i]=(alpha*pvec[i])+tmp; */
      /* } */

      scalar_xpy_d(rvec, -alpha, Av, rvec, N);

      /* tmp=0.0; */
      /* for(i=0;i<N;i++) */
      /* { */
      /*   tmp=rvec[i]; */
      /*   rvec[i]=(-alpha*Av[i])+tmp; */
      /* } */

      scalar_xpy_d(pvec, beta, pvec, rvec, N);

      /* tmp=0.0; */
      /* for(i=0;i<N;i++) */
      /* { */
      /*   tmp=rvec[i]; */
      /*   pvec[i]=(beta*pvec[i])+tmp; */
      /* } */

    }
  }

  free(Ar);
  free(Ap);
  free(delta);
  free(eta);
  free(zeta);
  free(rvec);
  free(pvec);
  free(Av);
  free(x_0);

  return 1;
}
