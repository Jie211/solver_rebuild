#include "gcr.h"

void GCR_init(double *rvec, double *Av, double *x_0, double *qq, double **qvec, double **pvec, const int N, const int restart)
{
  vec_init(rvec, 0.0, N);
  vec_init(Av, 0.0, N);
  vec_init(x_0, 0.0, N);
  vec_init(qq, 0.0,restart);
  vec_init_2(qvec, 0.0, N, restart);
  vec_init_2(pvec, 0.0, N, restart);
}

int GCR_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, const struct Parameter *para, const int N, const int NNZ, const bool f_isinner)
{
  int loop=0;
  int kloop, iloop;

  double *rvec, *Av, *x_0, *qq;
  double **qvec, **pvec;
  double alpha, beta, rnorm, bnorm;

  bool out_flag = false;
  double error = 0.0;
  double t_error = 0.0;

  int i_max = 0;
  double d_eps = 0.0;
  int i_restart=0;
  bool f_verbose=para->f_verbose;
  bool f_isVP=para->isVP;
  bool f_cuda=para->f_cuda;

  int exit_flag = 0.0;

  FILE *p_x = NULL, *p_his = NULL;

  double dot_tmp = 0.0;

  if(!f_isinner)
  {
    p_x = file_init("./output/GCR_x.txt", "w");
    p_his = file_init("./output/GCR_his.txt", "w");
  }

  if(!f_isinner)
  {
    i_restart = para->i_outer_restart;
  }else{
    i_restart = para->i_inner_restart;
  }

  if(f_cuda)
  {
    error_log("not done yet");
  }else{
    rvec = malloc_1d(N);
    Av = malloc_1d(N);
    x_0 = malloc_1d(N);
    qq = malloc_1d(i_restart);
    qvec = malloc_2d(N, i_restart);
    pvec = malloc_2d(N, i_restart);
  }

  //
  //start here
  //

  //init
  GCR_init(rvec, Av, xvec, qq, qvec, pvec, N, i_restart);

  if(f_isVP && f_isinner)
  {
    i_max = para->i_inner_maxloop;
    d_eps = para->d_inner_eps;
  }else{
    i_max = para->i_outer_maxloop;
    d_eps = para->d_outer_eps;
  }

  //backup xvc
  vec_copy(x_0, xvec, N);

  //b 2norm
  bnorm = norm_2_d(bvec, N);

  while(loop<i_max)
  {
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
    vec_copy(pvec[0], rvec, N);

    //Ap
    if(f_cuda)
    {
      error_log("not done yet");
    }else{
      MV_mult_CSR(qvec[0], val, col, ptr, pvec[0], N);
    }

    for(kloop=0;kloop<i_restart;kloop++)
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
#ifdef EBUG
        normal_log("GCR convergence");
#endif
        exit_flag = 1;
        out_flag=true;
        break;
      }else if(loop>=i_max)
      {
        exit_flag = 2;
        out_flag=true;
        break;
      }
      loop++;

      //(q, q)
      if(f_cuda)
      {
        error_log("not done");
        return -1;
      }else{
        dot_tmp = dot_d(qvec[kloop], qvec[kloop], N);
      }

      qq[kloop] = dot_tmp;


      //alpha = (r,q)/(q,q)
      if(f_cuda)
      {
        error_log("not done");
        return -1;
      }else{
        dot_tmp = dot_d(rvec, qvec[kloop], N);
      }

      alpha = dot_tmp / qq[kloop];

      //x=alpha*pvec[k]+xvec
      scalar_xpy_d(xvec, alpha, pvec[kloop], xvec, N);
      if(kloop == i_restart-1)
      {
        break;
      }

      //r=-alpha*qvec+rvec
      scalar_xpy_d(rvec, -alpha, qvec[kloop], rvec, N);

      //Ar
      if(f_cuda)
      {
        error_log("not done ");
        return -1;

      }else{
        MV_mult_CSR(Av, val, col, ptr, rvec, N);
      }

      //init p[k+1] q[k+1]
      vec_init(pvec[kloop+1], 0.0, N);
      vec_init(qvec[kloop+1], 0.0, N);

      for(iloop=0;iloop<=kloop;iloop++)
      {
        //beta=-(Az, qvec)/(q,q)
        if(f_cuda)
        {
          error_log("not done yet");
        }else{
          dot_tmp = dot_d(Av, qvec[iloop], N);
        }
        beta = -(dot_tmp)/qq[iloop];

        //pvec[k+1]=beta*pvec[i]+pvec[k+1]
        scalar_xpy_d(pvec[kloop+1], beta, pvec[iloop], pvec[kloop+1], N);
        //qvec[k+1]=beta*qvec[i]+qvec[k+1]
        scalar_xpy_d(qvec[kloop+1], beta, qvec[iloop], qvec[kloop+1], N);
      }
      //p[k]=z+p[k]
      vec_add_d(pvec[kloop+1], rvec, pvec[kloop+1], N);
      //q[k]=Az+q[k]
      vec_add_d(qvec[kloop+1], Av, qvec[kloop+1], N);
    }
    if(out_flag)
    {
      break;
    }
  }
  if(!f_isinner)
  {
    file_print(p_x, xvec, N);
    fclose(p_x);
    t_error = error_d_CRS(val, col, ptr, bvec, xvec, x_0, N);
    printf("|b-ax|2/|b|2=%.1f\n", t_error);
    printf("loop=%d\n", loop+1);
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
    free_1d(Av);
    free_1d(x_0);
    free_1d(qq);
    free_2d(qvec, i_restart);
    free_2d(pvec, i_restart);
  }
  if(f_isinner)
  {
    fclose(p_x);
    fclose(p_his);
  }
  return exit_flag;
}
