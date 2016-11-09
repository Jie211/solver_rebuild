#include "vpgmres.h"
void VPGMRES_Init(double *rvec, double *axvec, double *evec, double *vvec, double *vmtx, double *hmtx, double *yvec, double *wvec, double *avvec, double *hvvec, double *cvec, double *svec, double *x0vec, double *tmpvec, double *zmtx, double *zvec, double *xvec, const int N, const int i_restart)
{
  vec_init(rvec, 0.0, N);
  vec_init(axvec, 0.0, N);
  vec_init(vvec, 0.0, N);
  vec_init(wvec, 0.0, N);
  vec_init(avvec, 0.0, N);
  vec_init(x0vec, 0.0, N);
  vec_init(xvec, 0.0, N);
  vec_init(tmpvec, 0.0, N);
  vec_init(zvec, 0.0, N);
  
  vec_init(vmtx, 0.0, N*(i_restart+1));
  vec_init(hmtx, 0.0, N*(i_restart+1));
  vec_init(zmtx, 0.0, N*(i_restart+1));

  vec_init(yvec, 0.0, i_restart);
  vec_init(cvec, 0.0, i_restart);
  vec_init(svec, 0.0, i_restart);
  vec_init(evec, 0.0, i_restart);

  vec_init(hvvec, 0.0, i_restart*(i_restart+1));
}


int VPGMRES_CRS(double *val, int *col, int *ptr, double *Tval, int *Tcol, int *Tptr, double *bvec, double *xvec, struct Parameter *para, const int N, const int NNZ, const bool f_isinner){
  int i,j,k;
  FILE *p_x,*p_his;
  double *rvec, *axvec, *evec, *vvec, *vmtx, *hmtx, *yvec, *wvec, *avvec, *hvvec, *cvec, *svec, *x0vec, *tmpvec, *zmtx, *zvec, *x_0;
  double alpha,wv_ip;
  int count = 0;
  double tmp,tmp2,eps_now=0.0,b_norm;
  int flag = 0;
  int flag_break = 0;
  int error_message;
  double t_error;
  
  int ndata = N;

  int rs = para->i_outer_restart;
  int i_max = para->i_outer_maxloop;
  bool f_verbose = para->f_verbose;
  double d_eps = para->d_outer_eps;
  bool f_cuda = para->f_cuda;

  if(f_cuda)
  {
    error_log("not done yet");
  }else{
    rvec = malloc_1d(N);
    axvec = malloc_1d(N);
    evec = malloc_1d(rs);
    vvec = malloc_1d(N);
    vmtx = malloc_1d(N*(rs+1));
    hmtx = malloc_1d(N*(rs+1));
    yvec = malloc_1d(rs);
    wvec = malloc_1d(N);
    avvec = malloc_1d(N);
    hvvec = malloc_1d(rs*(rs+1));
    cvec = malloc_1d(rs);
    svec = malloc_1d(rs);
    x0vec = malloc_1d(N);
    tmpvec = malloc_1d(N);
    zmtx = malloc_1d(N*(rs+1));
    zvec = malloc_1d(N);
    x_0 = malloc_1d(N);
  }


  p_x = file_init("./output/VPGMRES_x.txt", "w");
  p_his = file_init("./output/VPGMRES_his.txt", "w");

  VPGMRES_Init(rvec, axvec, evec, vvec, vmtx, hmtx, yvec, wvec, avvec, hvvec, cvec, svec, x0vec, tmpvec, zmtx, zvec, xvec, N, rs);
  
  /* b_norm = vector_norm_2(bvec,ndata); */
  b_norm = norm_2_d(bvec, ndata);
  vec_copy(x_0, xvec, ndata);

  printf("@@@@@@@@@@MAX = %d\n", i_max);

  //outer loop 
  for(count=0;count<i_max;){
    //Ax0
    MV_mult_CSR(axvec, val, col, ptr, xvec, N);

    //r0 = b - Ax0
    vec_sub(rvec, bvec, axvec, N);

    //2norm(r)
    tmp = norm_2_d(rvec, ndata);

    //v0 = r0 / 2norm(r)
    //v0 >> vmtx[0][]
    for(i=0;i<ndata;i++){
      vvec[i] = rvec[i] / tmp;
      vmtx[0*ndata+i] = vvec[i];
    }


    evec[0] = tmp;

    for(i=1;i<rs;i++){
      evec[i]=0.0;
    }

    //inner loop k-> restart
    for(k=0;k<rs-1;k++){
      eps_now = fabs(evec[k]) / b_norm;
      if(f_verbose)
        printf("%d %.12e\n",count+1, eps_now);
      if(count >= i_max){
        flag_break = 1;
        break;
      }
      fprintf(p_his,"%d %.12e\n",count+1,eps_now);
      //if over eps break
      if(eps_now <= d_eps){
        solve_Hye(hmtx,yvec,evec,k,ndata);

        //epsilon yv
        vec_init(tmpvec, 0.0, ndata);

        for(i=0;i<k;i++){
          for(j=0;j<ndata;j++){
            tmpvec[j] += yvec[i] * zmtx[i*ndata+j];
          }
        }
        
        //x = x0 + epsilon yv
        /* for(i=0;i<ndata;i++){ */
        /*   xvec[i] = x0vec[i] + tmpvec[i]; */
        /* } */
        vec_add(xvec, x0vec, tmpvec, ndata);

        flag = 1;
        break;
      }
      //inner solver
      error_message = inner_selecter(para, vvec, zvec, val, col, ptr, Tval, Tcol, Tptr, N, NNZ);
      if(error_message == -1){
        error_log("error in vpgmres - inner_selecter");
      }

      /* //Av & W */
      MV_mult_CSR(wvec, val, col, ptr, zvec, N);

      for(i=0;i<ndata;i++){
        zmtx[k*ndata+i] = zvec[i];
      }

      //h_i_k & W  update
      for(i=0;i<=k;i++){
        wv_ip=0.0;
        for(j=0;j<ndata;j++){
          wv_ip+=wvec[j] * vmtx[i*ndata+j];
        }
        hmtx[i*ndata+k] = wv_ip;
        for(j=0;j<ndata;j++){
          wvec[j]=wvec[j] - wv_ip * vmtx[i*ndata+j];
        }
      }

      //h_k+1_k update
      tmp = norm_2_d(wvec, ndata);
      hmtx[(k+1)*ndata+k]=tmp;

      //v update
      for(i=0;i<ndata;i++){
        vvec[i] = wvec[i] / tmp;
        vmtx[(k+1)*ndata+i] = vvec[i];
      }

      //h_ update
      for(i=0;i<=(k-1);i++){
        tmp=hmtx[i*ndata+k];
        tmp2=hmtx[(i+1)*ndata+k];
        hmtx[i*ndata+k] = cvec[i] * tmp - svec[i] * tmp2;
        hmtx[(i+1)*ndata+k] = svec[i] * tmp + cvec[i] * tmp2;
      }

      //alpha = root(h_kk * h_kk + h_k+1_k * h_k+1_k)
      alpha = sqrt(hmtx[k*ndata+k] * hmtx[k*ndata+k] + hmtx[(k+1)*ndata+k] * hmtx[(k+1)*ndata+k]);

      cvec[k] = hmtx[k*ndata+k] / alpha;
      svec[k] = -hmtx[(k+1)*ndata+k] / alpha;
      evec[k+1] = svec[k] * evec[k];
      evec[k] = cvec[k] * evec[k];
      hmtx[k*ndata+k] = cvec[k] * hmtx[k*ndata+k] - svec[k] * hmtx[(k+1)*ndata+k];
      hmtx[(k+1)*ndata+k] = 0.0;

      count ++;
    }

    if(flag == 1){
      break;
    }
    if(flag_break == 1){
      break;
    }    
    solve_Hye(hmtx, yvec, evec, rs-1, ndata);

    for(i=0;i<ndata;i++){
      tmpvec[i]=0.0;
    }
    for(i=0;i<rs;i++){
      for(j=0;j<ndata;j++){
        tmpvec[j] += yvec[i] * zmtx[i*ndata+j];
      }
    }

    for(i=0;i<ndata;i++){
      xvec[i] = x0vec[i] + tmpvec[i];
    }

    for(i=0;i<ndata;i++){
      x0vec[i] = xvec[i];
    }
  }

  file_print(p_x, xvec, N);
  fclose(p_x);
  fclose(p_his);
  t_error = error_d_CRS(val, col, ptr, bvec, xvec, x_0, N);
  printf("|b-ax|2/|b|2=%.1f\n", t_error);
  printf("loop=%d\n", count+1);

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
    free_1d(zmtx);
    free_1d(zvec);
    free_1d(x_0);
  }
  if(flag){
    return 1;
  }else if(flag_break){
    return 2;
  }
  return -1;

}
