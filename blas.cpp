#include "blas.h"

double norm_1_d(double *v, const int N)
{
  int i;
  double tmp = 0.0;
  for(i=0;i<N;i++)
  {
    tmp += fabs(v[i]);
  }
  return tmp;
}

double norm_2_d(double *v, const int N)
{
  int i;
  double tmp = 0.0;
  for(i=0;i<N;i++)
  {
    tmp += v[i]*v[i];
  }
  return sqrt(tmp);
}

void MV_mult_CSR(double *out, const double *val, const int *col, const int *ptr, const double *vec, const int N)
{
  int i, j;
  double tmp = 0.0;
#pragma omp parallel for private(j) reduction(+:tmp) schedule(static) firstprivate(out, val, vec) lastprivate(out)
  for(i=0;i<N;i++){
    tmp=0.0;
    for(j=ptr[i];j<ptr[i+1];j++){
      tmp+=val[j] * vec[col[j]];
    }
    out[i]=tmp;
  }
}

void vec_sub(double *out, const double *x, const double *y, const int N)
{
  int i;
  for(i=0;i<N;i++)
  {
    out[i] = x[i] - y[i];
  }
}

void vec_add(double *out, const double *x, const double *y, const int N)
{
  int i;
  for(i=0;i<N;i++)
  {
    out[i] = x[i] + y[i];
  }
}

double dot_d(const double *x, const double *y, const int N)
{
  int i;
  double tmp = 0.0;
  for(i=0;i<N;i++)
  {
    tmp += x[i]*y[i];
  }
  return tmp;
}

void scalar_d(double *out, const double a, const double *x, const int N)
{
  int i;
  for(i=0;i<N;i++)
  {
    out[i] = a * x[i];
  }
}

void scalar_xpy_d(double *out, const double a, const double *x, const double *y, const int N)
{
  int i;
  double tmp;
  for(i=0;i<N;i++)
  {
    tmp = y[i];
    out[i] = (a * x[i]) + tmp;
  }
}

double error_d_CRS(double *val, const int *col, const int *ptr, const double *b, const double *x_new, const double *x_0, const int N)
{
  double tmp1=0.0;
  double tmp2=0.0;
  double tmp=0.0;
  double *Ax, *Ax_0;
  int i, j;

  Ax=(double *)malloc(sizeof(double)*N);
  Ax_0=(double *)malloc(sizeof(double)*N);

  for(i=0;i<N;i++)
  {
    tmp = 0.0;
    for(j=ptr[i];j<ptr[i+1];j++)
    {
      tmp += val[j] * x_new[col[j]];
    }
    Ax[i] = b[i] - tmp;
  }
  for(i=0;i<N;i++)
  {
    tmp = 0.0;
    for(j=ptr[i];j<ptr[i+1];j++)
    {
      tmp += val[j] * x_0[col[j]];
    }
    Ax_0[i] = b[i] - tmp;
  }
  tmp1 = norm_2_d(Ax, N);
  tmp2 = norm_2_d(Ax_0, N);
  tmp = log10(tmp1 / tmp2);

  free(Ax);
  free(Ax_0);

  return tmp;
}

void solve_Hye(double *h, double *y, double *e, const int n, const int size)
{
  int i, j;
  double tmp;
  for(i=0;i<n;i++)
  {
    y[i] = 0.0;
  }
  for(i=n-1;i>=0;i--)
  {
    tmp = 0.0;
    for(j=i+1;j<n;j++)
    {
      tmp += y[j] * h[i*size+j];
    }
    y[i] = (e[i] - tmp)/h[i*size+i];
  }
}

void cal_arap_kskipcg_d(double **Ar, double **Ap, double *val, int *col, int *ptr, double *rvec, double *pvec, const int N, const int kskip)
{
  int i, j, ii;
  double tmp1=0.0;
  double tmp2=0.0;
#pragma omp parallel for private(j) reduction(+:tmp1, tmp2) schedule(static) firstprivate(Ar, Ap, val, pvec, rvec) lastprivate(Ar, Ap)
  for(i=0;i<N;i++){
    tmp1=0.0;
    tmp2=0.0;
    for(j=ptr[i];j<ptr[i+1];j++){
      tmp1 += val[j]*rvec[col[j]];
      tmp2 += val[j]*pvec[col[j]];
    }
    Ar[0][i]=tmp1;
    Ap[0][i]=tmp2;
  }
  for(ii=1;ii<2*kskip+2;ii++){
#pragma omp parallel for private(i, j) reduction(+:tmp1, tmp2) schedule(static) firstprivate(Ar, Ap, val) lastprivate(Ar, Ap)
    for(i=0;i<N;i++){
      tmp1=0.0;
      tmp2=0.0;
      for(j=ptr[i];j<ptr[i+1];j++){
        if(ii<2*kskip){
          tmp1 += val[j]*Ar[(ii-1)][col[j]];
        }
        tmp2 += val[j]*Ap[(ii-1)][col[j]];
      }
      if(ii<2*kskip){
        Ar[ii][i]=tmp1;
      }
      Ap[ii][i]=tmp2;
    }
  }
}

void cal_arap_kskipbicg_d(double **Ar, double **Ap, double *val, int *col, int *ptr, double *rvec, double *pvec, const int N, const int kskip)
{
  int i, j, ii;
  double tmp1=0.0;
  double tmp2=0.0;
#pragma omp parallel for private(j) reduction(+:tmp1, tmp2) schedule(static) firstprivate(Ar, Ap, val, pvec, rvec) lastprivate(Ar, Ap)
  for(i=0;i<N;i++){
    tmp1=0.0;
    tmp2=0.0;
    for(j=ptr[i];j<ptr[i+1];j++){
      tmp1 += val[j]*rvec[col[j]];
      tmp2 += val[j]*pvec[col[j]];
    }
    Ar[0][i]=tmp1;
    Ap[0][i]=tmp2;
  }
  for(ii=1;ii<2*kskip+2;ii++){
#pragma omp parallel for private(i, j) reduction(+:tmp1, tmp2) schedule(static) firstprivate(Ar, Ap, val) lastprivate(Ar, Ap)
    for(i=0;i<N;i++){
      tmp1=0.0;
      tmp2=0.0;
      for(j=ptr[i];j<ptr[i+1];j++){
        if(ii<2*kskip+1){
          tmp1 += val[j]*Ar[(ii-1)][col[j]];
        }
        tmp2 += val[j]*Ap[(ii-1)][col[j]];
      }
      if(ii<2*kskip+1){
        Ar[ii][i]=tmp1;
      }
      Ap[ii][i]=tmp2;
    }
  }
}


void cal_deltaetazeta_kskipcg_d(double *delta, double *eta, double *zeta, double **Ar, double **Ap, double *rvec, double *pvec, const int N, const int kskip)
{
  int i, j;
  double tmp1=0.0;
  double tmp2=0.0;
  double tmp3=0.0;
#pragma omp parallel for private(j) reduction(+:tmp1, tmp2, tmp3) schedule(static) firstprivate(delta, eta, zeta, Ar, rvec, Ap, pvec) lastprivate(delta, eta, zeta)
  for(i=0;i<2*kskip+2;i++){
    tmp1=0.0;
    tmp2=0.0;
    tmp3=0.0;
    for(j=0;j<N;j++){
      if(i<2*kskip){
        tmp1 += rvec[j]*Ar[i][j];
      }
      if(i<2*kskip+1){
        tmp2 += rvec[j]*Ap[i][j];
      }
      tmp3 += pvec[j]*Ap[i][j];
    }
    if(i<2*kskip){
      delta[i]=tmp1;
    }
    if(i<2*kskip+1){
      eta[i]=tmp2;
    }
    zeta[i]=tmp3;
  }
}

void cal_theta_eta_rho_phi_kskipcg_d(double *theta, double *eta, double *rho, double *phi, double **Ar, double **Ap, double *rvec, double *pvec, double *r_vec, double *p_vec, const int N, const int kskip)
{
  int i, j;
    /* theta (2*i_kskip); */
    /* eta = (2*i_kskip+1); */
    /* rho = (2*i_kskip+1); */
    /* phi = (2*i_kskip+2); */
    //theta = (r*, Ar)
    //eta = (r*, Ap)
    //rho = (p*, Ar)
    //phi = (p*, Ap)
  double tmp1=0.0;//theta
  double tmp2=0.0;//eta
  double tmp3=0.0;//rho
  double tmp4=0.0;//phi
#pragma omp parallel for private(j) reduction(+:tmp1, tmp2, tmp3, tmp4) schedule(static) firstprivate(theta, eta, rho, phi, Ar, rvec, Ap, pvec, r_vec, p_vec) lastprivate(theta, eta, rho, phi)
  for(i=0;i<2*kskip+2;i++){
    tmp1=0.0;
    tmp2=0.0;
    tmp3=0.0;
    tmp4=0.0;
    for(j=0;j<N;j++){
      if(i<2*kskip){
        tmp1 += r_vec[j]*Ar[i][j];
      }
      if(i<2*kskip+1){
        tmp2 += r_vec[j]*Ap[i][j];
        tmp3 += p_vec[j]*Ar[i][j];
      }
      tmp4 += p_vec[j]*Ap[i][j];
    }
    if(i<2*kskip){
      theta[i]=tmp1;
    }
    if(i<2*kskip+1){
      eta[i]=tmp2;
      rho[i]=tmp3;
    }
    phi[i]=tmp4;
  }
}

void cal_arap_kskipcr_d(double *Ar, double *Ap, double *val, int *col, int *ptr, double *rvec, double *pvec, const int N, const int kskip)
{
  int i, j, ii;
  double tmp1=0.0;
  double tmp2=0.0;

#pragma omp parallel for private(j) reduction(+:tmp1, tmp2) schedule(static) firstprivate(Ar, Ap, val, pvec, rvec) lastprivate(Ar, Ap)
  for(i=0;i<N;i++)
  {
    tmp1=0.0;
    tmp2=0.0;
    for(j=ptr[i];j<ptr[i+1];j++)
    {
      tmp1 += val[j] * rvec[col[j]];
      tmp2 += val[j] * pvec[col[j]];
    }
    Ar[i]=tmp1;
    Ap[i]=tmp2;
  }

#pragma omp parallel for private(i, j) reduction(+:tmp1, tmp2) schedule(static) firstprivate(Ar, Ap, val) lastprivate(Ar, Ap)
  for(ii=1;ii<2*kskip+1;ii++){
    for(i=0;i<N;i++){
      tmp1=0.0;
      tmp2=0.0;
      for(j=ptr[i];j<ptr[i+1];j++){
        tmp1 += val[j]*Ar[(ii-1)*N+col[j]];
        tmp2 += val[j]*Ap[(ii-1)*N+col[j]];
      }
      Ar[ii*N+i]=tmp1;
      Ap[ii*N+i]=tmp2;
    }
  }
}

void cal_deltaetazeta_kskipcr_d(double *delta, double *eta, double *zeta, double *Ar, double *Ap, double *rvec, const int N, const int kskip)
{
  int i, j;
  double tmp1=0.0;
  double tmp2=0.0;
  double tmp3=0.0;
#pragma omp parallel for private(j) reduction(+:tmp1, tmp2, tmp3) schedule(static) firstprivate(delta, eta, zeta, Ar, rvec, Ap) lastprivate(delta, eta, zeta)
  for(i=0;i<2*kskip+1;i++){
    tmp1=0.0;
    tmp2=0.0;
    tmp3=0.0;
    for(j=0;j<N;j++){
      tmp1 += rvec[j]*Ar[i*N+j];
      tmp2 += Ap[0*N+j]*Ap[i*N+j];
      tmp3 += rvec[j]*Ap[i*N+j];
    }
    delta[i]=tmp1;
    eta[i]=tmp2;
    zeta[i]=tmp3;
  }
}

void Transpose_d(double *val, int *col, int *ptr, double *Tval, int *Tcol, int *Tptr, const int N, const int NNZ){

  int col_counter = 0;
  memset(Tptr, -1, sizeof(int)*(N+1));

  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      for(int k=ptr[j];k<ptr[j+1];k++){
        if(col[k]==i){
          if(Tptr[i]==-1){
            Tptr[i] = col_counter;
          }
          Tcol[col_counter] = j;
          Tval[col_counter] = val[k];
          col_counter++;
          continue;
        }
      }
    }
  }

  Tptr[N] = NNZ;
}
