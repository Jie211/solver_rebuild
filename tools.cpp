#include "tools.h"

void error_log(char *output)
{
  printf("\x1b[31m");
  printf("[×] %s\n",output);
  printf("\x1b[0m");
}

void warning_log(char *output)
{
  printf("\x1b[33m");
  printf("[△] %s\n",output);
  printf("\x1b[0m");
}

void normal_log(char *output)
{
  printf("\x1b[32m");
  printf("[○] %s\n",output);
  printf("\x1b[0m");
}


double *malloc_1d(int num)
{
  double *tmp;
  tmp=(double *)malloc(sizeof(double)*num);
  if(!tmp)
  {
    error_log((char*)"1 dim double malloc error");
    exit(-1);
  }
  return tmp;
}

double **malloc_2d(int num_x, int num_y)
{
  double **tmp;
  int i;
  tmp=(double **)malloc(sizeof(double*)*num_y);
  if(!tmp)
  {
    error_log((char*)"2 dim double malloc error[1]");
    exit(-1);
  }
  for(i=0; i<num_y; i++)
  {
    tmp[i] = (double *)malloc(sizeof(double)*num_x);
    if(!tmp[i])
    {
      error_log((char*)"2 dim double malloc error[2]");
      exit(-1);
    }
  }
  return (tmp);
}

int *malloc_1i(int num)
{
 int *tmp;
 tmp=(int *)malloc(sizeof(int)*num);
  if(!tmp)
  {
    error_log((char*)"1 dim int malloc error");
    exit(-1);
  }
  return tmp;
}

int** malloc_2i(int num_x, int num_y)
{
  int **tmp, i;
  tmp=(int **)malloc(sizeof(int*)*num_y);
  if(!tmp)
  {
    error_log((char*)"2 dim double malloc error[1]");
    exit(-1);
  }
  for(i=0; i<num_y; i++)
  {
    tmp[i] = (int *)malloc(sizeof(int)*num_x);
    if(!tmp[i])
    {
      error_log((char*)"2 dim double malloc error[2]");
      exit(-1);
    }
  }
  return tmp;
}

void free_1d(double *ptr)
{
  free(ptr);
}

void free_2d(double **ptr, int num)
{
  int i;
  for(i=0;i<num;i++)
  {
    free(ptr[i]);
  }
  free(ptr);
}

void free_1i(int *ptr)
{
  free(ptr);
}

void free_2i(int **ptr, int num)
{
  int i;
  for(i=0;i<num;i++)
  {
    free(ptr[i]);
  }
  free(ptr);
}

int error_handle(int error_code, char *msg)
{
  if(error_code!=0)
  {
    error_log(msg);
    return -1;
  }else{
    return 0;
  }
}

int set_openmp_thread(const int thread)
{
  int error;
  char env[128];
  char setThreads[10];
  sprintf(setThreads, "%d", thread);
  error = setenv("OMP_NUM_THREADS", setThreads, 1);
  if(error!=0)
    return -1;
  strcpy(env, "OMP_NUM_THREADS=");
  strcat(env, setThreads);
  error = putenv(env);
  if(error!=0)
    return -1;
  omp_set_num_threads(thread);
#ifdef EBUG
  normal_log((char*)"pass set openmp thread");
#endif
  return 0;
}

FILE *file_init(const char *name, const char *mode)
{
  FILE *tmp;
  if((tmp = fopen(name, mode)) == NULL)
  {
    warning_log((char*)"file open failed");
    exit(-1);
  }
  return (tmp);
}

void vec_init(double *v, const double val, const int size)
{
  int i;
  for(i=0;i<size;i++)
  {
    v[i]=val;
  }
}

void vec_init_2(double **v, const double val, const int sizex, const int sizey)
{
  int i, j;
  for(i=0;i<sizey;i++)
  {
    for(j=0;j<sizex;j++)
    {
      v[i][j] = val;
    }
  }
}

void vec_copy(double *t, double *f, const int size)
{
  int i;
  for(i=0;i<size;i++)
  {
    t[i] = f[i];
  }
}

void file_print(FILE *fp, double *vec, const int N)
{
  int i;
  for(i=0;i<N;i++)
  {
    fprintf(fp, "%d %.12e\n", i, vec[i]);
  }
}

void vec_add_d(double *out, double *x, double *y, const int N)
{
  int i;
  for(i=0;i<N;i++)
  {
    out[i] = x[i] + y[i];
  }
}
