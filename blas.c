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
