#include "kskipcg.h"

void KSKIPBICG_init(double **Ap, double **Ar, double *theta, double *eta, double *rho, double *phi, double *rvec, double *pvec, double *r_vec, double *p_vec, double *Av, double *xvec, const int N, const int i_kskip)
{
  vec_init_2(Ar, 0.0, N, (2*i_kskip+1));
  vec_init_2(Ap, 0.0, N, (2*i_kskip+2));
  vec_init(theta, 0.0, 2*i_kskip);
  vec_init(eta, 0.0, 2*i_kskip+1);
  vec_init(rho, 0.0, 2*i_kskip+1);
  vec_init(phi, 0.0, 2*i_kskip+2);
  vec_init(rvec, 0.0, N);
  vec_init(pvec, 0.0, N);
  vec_init(r_vec, 0.0, N);
  vec_init(p_vec, 0.0, N);
  vec_init(Av, 0.0, N);
  vec_init(xvec, 0.0, N);
}
int KSKIPBICG_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, const struct Parameter *para, const int N, const int NNZ, const bool f_isinner)
{
  int nloop, iloop, jloop;

  double *theta, *eta, *rho, *phi;
  double *rvec, *pvec, *r_vec, *p_vec, *Av, *x_0;
  double **Ap, **Ar;

  double rnorm, bnorm, alpha, beta, gamma;

  int exit_flag = 2;
  double error = 0.0;

  double t_error;

  FILE *p_x=NULL, *p_his=NULL;

  int i_max = 0;
  double d_eps = 0.0;
  int i_kskip = 0;
  bool f_cuda = para->f_cuda;
  /* int i_fix; */
  bool f_verbose = para->f_verbose;

  if(!f_isinner)
  {
    p_x=file_init("./output/KSKIPBICG_x.txt", "w");
    p_his=file_init("./output/KSKIPBICG_his.txt", "w");
  }

  if(!f_isinner)
  {
    i_max = para->i_outer_maxloop;
    i_kskip = para->i_outer_kskip;
    /* i_fix = para->i_outer_fix; */
    d_eps = para->d_outer_eps;
  }else{
    i_max = para->i_inner_maxloop;
    i_kskip = para->i_inner_kskip;
    /* i_fix = para->i_inner_fix; */
    d_eps = para->d_inner_eps;
  }

  if(f_cuda)
  {
    error_log("not done yet");
  }else{
    Ar=malloc_2d(N, 2*i_kskip+1);
    Ap=malloc_2d(N, 2*i_kskip+2);
    
    theta = malloc_1d(2*i_kskip);
    eta = malloc_1d(2*i_kskip+1);
    rho = malloc_1d(2*i_kskip+1);
    phi = malloc_1d(2*i_kskip+2);
    

    rvec=malloc_1d(N);
    pvec=malloc_1d(N);
    r_vec=malloc_1d(N);
    p_vec=malloc_1d(N);

    Av=malloc_1d(N);
    x_0=malloc_1d(N);
  }
 

  int i, j, k;
  double *Tval;
  int *Tcol, *Tptr;
  int Tlock=0;
  int col_counter=0;
  Tval=(double *)malloc(sizeof(double)*NNZ);
  Tcol=(int *)malloc(sizeof(int)*NNZ);
  Tptr=(int *)malloc(sizeof(int)*N+1);

  for(i=0;i<N;i++){
    Tlock=0;
    for(j=0;j<N;j++){
      for(k=ptr[j];k<ptr[j+1];k++){
        if(col[k]==i){
          if(Tlock==0){
            Tptr[i] = col_counter;
            Tlock=1;
          }
          Tcol[col_counter] = j;
          Tval[col_counter] = val[k];
          col_counter++;
          continue;
        }
      }
    }
  }
  Tptr[N]=NNZ;



  //
  //start here
  //

  //init
  KSKIPBICG_init(Ap, Ar, theta, eta, rho, phi, rvec, pvec, r_vec, p_vec, Av, xvec, N, i_kskip);

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

  //r*=r
  vec_copy(r_vec, rvec, N);

  //p*=r*
  vec_copy(p_vec, r_vec, N);

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

    //Ar-> Ar^2k+1
    //Ap-> Ap^2k+2
    if(f_cuda)
    {
      error_log("not done yet 2");
    }else{
      cal_arap_kskipbicg_d(Ar, Ap, val, col, ptr, rvec, pvec, N, i_kskip);
    }

    //gamma=(r*,r)
    if(f_cuda)
    {
      error_log("not done yet 3");
    }else{
      gamma = dot_d(r_vec, rvec, N);
    }

    /* theta (2*i_kskip); */
    /* eta = (2*i_kskip+1); */
    /* rho = (2*i_kskip+1); */
    /* phi = (2*i_kskip+2); */
    //theta = (r*, Ar)
    //eta = (r*, Ap)
    //rho = (p*, Ar)
    //phi = (p*, Ap)
    if(f_cuda)
    {
      error_log("not done yet 4");
    }else{
      cal_theta_eta_rho_phi_kskipcg_d(theta, eta, rho, phi, Ar, Ap, rvec, pvec, r_vec, p_vec, N, i_kskip);
    }

    for(iloop=nloop;iloop<=nloop+i_kskip;iloop++)
    {
      //alpha = gamma/phi_1
      alpha=gamma/phi[0];
      //beta = 1 - (eta_1 + rho_1 - alpha*phi_2)/phi_1
      beta=1.0 - (eta[0] + rho[0] - alpha*phi[1])/phi[0];
      //gamma = gamma - alpha*eta_1 - alpha*rho_1 + alpha^2*phi_2
      /* gamma = gamma - alpha*eta[0] - alpha*rho[0] + alpha*alpha*phi[1]; */
      gamma = gamma*beta;

      //update theta eta rho phi
      for(jloop=0; jloop<2*i_kskip-2*(iloop-nloop);jloop++){
        theta[jloop] = theta[jloop] - alpha*eta[jloop+1] - alpha*rho[jloop+1] + alpha*alpha*phi[jloop+2];
        double tmp = theta[jloop] - alpha*beta*phi[jloop+1];
        eta[jloop] = tmp + beta*eta[jloop];
        rho[jloop] = tmp + beta*rho[jloop];
        phi[jloop] = eta[jloop] + rho[jloop] - theta[jloop] + beta*beta*phi[jloop];
      }

      //Ap
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

      //A^Tp*
      if(f_cuda)
      {
        error_log("not done yet 5");
      }else{
        MV_mult_CSR(Av, Tval, Tcol, Tptr, p_vec, N);
      }

      //r*=r*-alpha*A^Tp*
      scalar_xpy_d(r_vec, -alpha, Av, r_vec, N);

      //p=r+beta*p
      scalar_xpy_d(pvec, beta, pvec, rvec, N);

      //p*=r*+beta*p*
      scalar_xpy_d(p_vec, beta, p_vec, r_vec, N);
    }
  }
  if(!f_isinner)
  {
    file_print(p_x, xvec, N);
    fclose(p_x);
    fclose(p_his);
    t_error = error_d_CRS(val, col, ptr, bvec, xvec, x_0, N);
    printf("|b-ax|2/|b|2=%.1f\n", t_error);
    printf("loop=%d\n", nloop+1);
  }else{
    /* fclose(p_x); */
    /* fclose(p_his); */
  }
  if(f_isinner && f_verbose)
  {
    printf("Inner %d %.12e\n", nloop+1, error);
  }
  if(f_cuda)
  {
    error_log("not done yet 6");
  }else{
    free_1d(theta);
    free_1d(eta);
    free_1d(rho);
    free_1d(phi);
    free_1d(rvec);
    free_1d(pvec);
    free_1d(r_vec);
    free_1d(p_vec);
    free_1d(Av);
    free_1d(x_0);
    free_2d(Ap, 2*i_kskip+2);
    free_2d(Ar, 2*i_kskip+1);
    free(Tval);
    free(Tcol);
    free(Tptr);
  }
  
  return exit_flag;
}
