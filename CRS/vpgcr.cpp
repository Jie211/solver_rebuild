#include "vpgcr.h"

void VPGCR_Init(double *rvec, double *zvec, double *Av, double **qvec, double **pvec, double *qq, double *xvec, const int N, const int restart)
{
  vec_init(rvec, 0.0, N);
  vec_init(zvec, 0.0, N);
  vec_init(Av, 0.0, N);
  vec_init(qq, 0.0, N);
  vec_init(xvec, 0.0, N);
  vec_init_2(qvec, 0.0, N, restart);
  vec_init_2(pvec, 0.0, N, restart);
}

int VPGCR_CRS(double *val, int *col, int *ptr, double *Tval, int *Tcol, int *Tptr, double *bvec, double *xvec, struct Parameter *para, const int N, const int NNZ, const bool f_isinner)
{
  int loop=0, kloop, iloop;

  double *rvec, *zvec, *Av;
  double **qvec, **pvec, *qq, *x_0;
  double alpha, beta, rnorm, bnorm;
  double dot_tmp;

  int i_max = para->i_outer_maxloop;
  double d_eps = para->d_outer_eps;
  bool f_cuda = para->f_cuda;
  int i_restart = para->i_outer_restart;
  bool f_verbose = para->f_verbose;

  int exit_flag = 2;
  bool out_flag = false;
  int get_error = 0.0;
  double t_error = 0.0;
  double error = 0.0;
  
  FILE *p_x = NULL, *p_his = NULL;

  p_x = file_init("./output/VPGCR_x.txt", "w");
  p_his = file_init("./output/VPGCR_his.txt", "w");

  if(f_cuda)
  {
    error_log((char*)"not done yet");
  }else{
    rvec = malloc_1d(N);
    zvec = malloc_1d(N);
    Av = malloc_1d(N);
    qq = malloc_1d(N);
    x_0 = malloc_1d(N);
    qvec = malloc_2d(N, i_restart);
    pvec = malloc_2d(N, i_restart);
  }

  VPGCR_Init(rvec, zvec, Av, qvec, pvec, qq, xvec, N, i_restart);

  //
  //start
  //

  //backup xvec
  vec_copy(x_0, xvec, N);

  //b 2norm
  bnorm = norm_2_d(bvec, N);

  while(loop<i_max){
    //Ax
    if(f_cuda)
    {
      error_log((char*)"not done yet");
    }else{
      MV_mult_CSR(Av, val, col, ptr, xvec, N);
    }

    //r=b-Ax
    vec_sub(rvec, bvec, Av, N);

    //init pvec[0]
    vec_init(pvec[0], 0.0, N);

    //solve p by Ap[0]=r
    get_error = inner_selecter(para, rvec, pvec[0], val, col, ptr, Tval, Tcol, Tptr, N, NNZ);
    if(get_error == -1)
    {
      error_log((char*)"error in vpgcr - inner_selecter");
    }

    //q[0*ndata+x]=A*p[0*ndata+x]
    if(f_cuda)
    {
      error_log((char*)"cuda not done yet");
    }else{
      MV_mult_CSR(qvec[0], val, col, ptr, pvec[0], N);
    }

    for(kloop=0;kloop<i_restart;kloop++){

      rnorm = norm_2_d(rvec, N);

      error=rnorm/bnorm;

      if(f_verbose){
        printf("Outer %d %.12e\n", loop+1, error);
      }
      fprintf(p_his, "%d %.12e\n", loop+1, error);
      if(error <= d_eps)
      {
        exit_flag = 1;
        out_flag = true;
        break;
      }else if(loop >= i_max){
        exit_flag = 2;
        out_flag = true;
        break;
      }
      loop++;

      //(q,q)
      if(f_cuda)
      {
        error_log((char*)"not done yet");
      }else{
        qq[kloop] = dot_d(qvec[kloop], qvec[kloop], N);
      }
      

      //alpha=(r, q)/(q, q)
      if(f_cuda)
      {
        error_log((char*)"not done yet");
      }else{
        dot_tmp = dot_d(rvec, qvec[kloop], N);
      }
      alpha = dot_tmp / qq[kloop];

      //x=alpha*pvec[k]+xvec
      scalar_xpy_d(xvec, alpha, pvec[kloop], xvec, N);

      if(kloop==i_restart-1){
        break;
      }

      //r=-alpha*qvec+rvec
      scalar_xpy_d(rvec, -alpha, qvec[kloop], rvec, N);

      //init z
      vec_init(zvec, 0.0, N);

      //Az=r
      get_error = inner_selecter(para, rvec, zvec, val, col, ptr, Tval, Tcol, Tptr, N, NNZ);
      if(get_error == -1)
      {
        error_log((char*)"error in vpgcr - inner_selecter");
      }


      //Az
      MV_mult_CSR(Av, val, col, ptr, zvec, N);

      //init p[k+1] q[k+1]
      vec_init(pvec[kloop+1], 0.0, N);
      vec_init(qvec[kloop+1], 0.0, N);

      for(iloop=0; iloop<=kloop; iloop++){
        //beta=-(Az,qvec)/(q,q)
        if(f_cuda){
          error_log((char*)"not done yet");
        }else{
          dot_tmp = dot_d(Av, qvec[iloop], N);
        }
        beta = -(dot_tmp / qq[iloop]);

        //pvec[k+1]=beta*pvec[i]+pvec[k+1]
        scalar_xpy_d(pvec[kloop+1], beta, pvec[iloop], pvec[kloop+1], N);

        //qvec[k+1]=beta*qvec[i]+qvec[k+1]
        scalar_xpy_d(qvec[kloop+1], beta, qvec[iloop], qvec[kloop+1], N);
      }

      //p[k]=z+p[k]
      vec_add(pvec[kloop+1], zvec, pvec[kloop+1], N);

      //q[k]=Az+q[k]
      vec_add(qvec[kloop+1], Av, qvec[kloop+1], N);
    }
    if(out_flag){
      break;
    }
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
    free_1d(zvec);
    free_1d(Av);
    free_1d(x_0);
    free_1d(qq);

    free_2d(qvec, i_restart);
    free_2d(pvec, i_restart);
  }
  return exit_flag;
}
