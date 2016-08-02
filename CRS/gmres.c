#include "gmres.h"

void GMRES_init(double *rvec, double *axvec, double *evec, double *vvec, double *vmtx, double *hmtx, double *yvec, double *wvec, double *avvec, double *hvvec, double *cvec, double *svec, double *x0vec, double *tmpvec, double *x_0, double *xvec, const int N, const int restart)
{
  int i;
  for(i=0;i<N;i++)
  {
    rvec[i]=0.0;
    axvec[i]=0.0;
    vvec[i]=0.0;
    wvec[i]=0.0;
    avvec[i]=0.0;
    x0vec[i]=0.0;
    xvec[i]=0.0;
    tmpvec[i]=0.0;
    x_0[i]=0.0;
  }
  for(i=0;i<N*(restart+1);i++)
  {
    vmtx[i]=0.0;
    hmtx[i]=0.0;
  }
  for(i=0;i<restart;i++)
  {
    yvec[i]=0.0;
    cvec[i]=0.0;
    svec[i]=0.0;
    evec[i]=0.0;
  }
  for(i=0;i<restart*(restart+1);i++)
  {
    hvvec[i]=0.0;
  }
}

int GMRES_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, const struct Parameter *para, const int N, const int NNZ, const bool f_isinner)
{
  int i, j, k;
  double *rvec, *axvec, *evec, *vvec, *vmtx, *hmtx, *yvec, *wvec, *avvec, *hvvec, *cvec, *svec, *x0vec, *tmpvec, *x_0;
  double alpha, wv_ip;
  int count = 0;
  double tmp, tmp2, error = 0.0, b_norm;

  int exit_flag = 2;
  double t_error;

  int i_restart=0;
  int i_max=0;
  double d_eps=0.0;
  /* bool f_isVP=para->isVP; */
  bool f_verbose=para->f_verbose;
  bool f_cuda=para->f_cuda;

  FILE *p_x=NULL, *p_his=NULL;

  if(!f_isinner)
  {
    p_x=file_init("./output/GMRES_x.txt", "w");
    p_his=file_init("./output/GMRES_his.txt", "w");
  }

  if(!f_isinner)
  {
    i_restart = para->i_outer_restart;
    i_max = para->i_outer_maxloop;
    d_eps = para->d_outer_eps;
  }else{
    i_restart = para->i_inner_restart;
    i_max = para->i_inner_maxloop;
    d_eps = para->d_inner_eps;
  }

  if(f_cuda)
  {
    error_log("not done yet");
  }else{
    rvec = malloc_1d(N);
    axvec= malloc_1d(N);
    evec = malloc_1d(i_restart);
    vvec = malloc_1d(N);
    vmtx = malloc_1d(N*(i_restart+1));
    hmtx = malloc_1d(N*(i_restart+1));
    yvec = malloc_1d(i_restart);
    wvec = malloc_1d(N);
    avvec= malloc_1d(N);
    hvvec= malloc_1d(i_restart*(i_restart+1));
    cvec = malloc_1d(i_restart);
    svec = malloc_1d(i_restart);
    x0vec= malloc_1d(N);
    tmpvec=malloc_1d(N);
    x_0  = malloc_1d(N);
  }

  //
  //start here
  //

  //init
  GMRES_init(rvec, axvec, evec, vvec, vmtx, hmtx, yvec, wvec, avvec, hvvec, cvec, svec, x0vec, tmpvec, x_0, xvec, N, i_restart);

  //b norm
  b_norm = norm_2_d(bvec, N);

  //backup xvec
  vec_copy(x_0, xvec, N);

  for(count=0;count<i_max;)
  {
    //Ax0
    if(f_cuda)
    {
      error_log("not done yet");
    }else{
      MV_mult_CSR(axvec, val, col, ptr, xvec, N);
    }

    //r0 = b-Ax0
    vec_sub(rvec, bvec, axvec, N);

    //2norm rvec
    tmp = norm_2_d(rvec, N);

    for(i=0;i<N;i++)
    {
      vvec[i] = rvec[i] / tmp;
      vmtx[0*N+i] = vvec[i];
    }

    vec_init(evec, 0.0, i_restart);
    evec[0] = tmp;

    for(k=0;k<i_restart-1;k++)
    {
      error = fabs(evec[k]) / b_norm;
      if(!f_isinner)
      {
        if(f_verbose)
        {
          printf("%d %.12e\n", count+1, error);
        }
        fprintf(p_his, "%d %.12e\n", count+1, error);
      }
      if(error <= d_eps)
      {
        solve_Hye(hmtx, yvec, evec, k, N);

        vec_init(tmpvec, 0.0, N);

        for(i=0;i<k;i++)
        {
          for(j=0;j<N;j++)
          {
            tmpvec[j] += yvec[i] * vmtx[i*N+j];
          }
        }

        vec_add(xvec, x0vec, tmpvec, N);

#ifdef EBUG
        normal_log("GCR convergence");
#endif
        exit_flag = 1;
        break;
      }

      //Av & w
      for(i=0;i<N;i++)
      {
        tmp = 0.0;
        for(j=ptr[i];j<ptr[i+1];j++)
        {
          tmp += val[j] * vmtx[k*N+col[j]];
        }
        avvec[i] = tmp;
        wvec[i] = avvec[i];
      }

      //h_i_k & W update
      for(i=0;i<=k;i++)
      {
        for(j=0;j<N;j++)
        {
          tmpvec[j] = vmtx[i*N+j];
        }
        wv_ip = 0.0;
        for(j=0;j<N;j++)
        {
          wv_ip+=wvec[j] * vmtx[i*N+j];
        }
        hmtx[i*N+k] = wv_ip;
        for(j=0;j<N;j++)
        {
          wvec[j] = wvec[j] - wv_ip * vmtx[i*N+j];
        }
      }

      //h_k+1 update
      tmp = norm_2_d(wvec, N);
      hmtx[(k+1)*N+k]=tmp;

      //v update
      for(i=0;i<N;i++)
      {
        vvec[i] = wvec[i] / tmp;
        vmtx[(k+1)*N+i] = vvec[i];
      }

      //h update
      for(i=0;i<=(k-1);i++)
      {
        tmp = hmtx[i*N+k];
        tmp2= hmtx[(i+1)*N+k];
        hmtx[i*N+k] = cvec[i] * tmp - svec[i] * tmp2;
        hmtx[(i+1)*N+k] = svec[i] * tmp + cvec[i] * tmp2;
      }

      //alpha = root(h_kk * h_kk + h_k+1_k * h_k+1_k)
      alpha = sqrt(hmtx[k*N+k] * hmtx[k*N+k] + hmtx[(k+1)*N+k] * hmtx[(k+1)*N+k]);

      cvec[k] = hmtx[k*N+k] / alpha;
      svec[k] = -hmtx[(k+1)*N+k] / alpha;
      evec[k+1] = svec[k] * evec[k];
      evec[k] = cvec[k] * evec[k];
      hmtx[k*N+k] = cvec[k] * hmtx[k*N+k] - svec[k] * hmtx[(k+1)*N+k];
      hmtx[(k+1)*N+k] = 0.0;

      count++;
    }

    if(exit_flag == 1)
    {
      break;
    }

    solve_Hye(hmtx, yvec, evec, i_restart-1, N);

    vec_init(tmpvec, 0.0, N);

    for(i=0;i<i_restart;i++)
    {
      for(j=0;j<N;j++)
      {
        tmpvec[j] += yvec[i] * vmtx[i*N+j];
      }
    }

    vec_add(xvec, x0vec, tmpvec, N);
    
    for(i=0;i<N;i++)
    {
      x0vec[i]=xvec[i];
    }
  }
  if(!f_isinner)
  {
    file_print(p_x, xvec, N);
    fclose(p_x);
    t_error = error_d_CRS(val, col, ptr, bvec, xvec, x_0, N);
    printf("|b-ax|2/|b|2=%.1f\n", t_error);
    printf("loop=%d\n", count+1);
  }
  if(f_isinner && f_verbose)
  {
    printf("Inner %d %.12e\n", count+1, error);
  }
  if(f_cuda)
  {
    error_log("not done yet");
  }else{
    free_1d(rvec);
    free_1d(axvec);
    free_1d(evec);
    free_1d(vvec);
    free_1d(vmtx);
    free_1d(hmtx);
    free_1d(yvec);
    free_1d(wvec);
    free_1d(avvec);
    free_1d(hvvec);
    free_1d(cvec);
    free_1d(svec);
    free_1d(x0vec);
    free_1d(tmpvec);
    free_1d(x_0);
  }
  if(f_isinner)
  {
    fclose(p_x);
    fclose(p_his);
  }
  return exit_flag;
}
