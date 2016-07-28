#include "tools.h"

void error_log(char *output)
{
  printf("\x1b[31m");
  printf("[×] %s\n",output);
}

void warning_log(char *output)
{
  printf("\x1b[33m");
  printf("[△] %s\n",output);
}

void normal_log(char *output)
{
  printf("\x1b[32m");
  printf("[◯] %s\n",output);
}

double *malloc_1d(int num)
{
  double *tmp;
  tmp=(double *)malloc(sizeof(double)*num);
  if(!tmp)
  {
    error_log("1 dim double malloc error");
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
    error_log("2 dim double malloc error[1]");
    exit(-1);
  }
  for(i=0; i<num_y; i++)
  {
    tmp[i] = (double *)malloc(sizeof(double)*num_x);
    if(!tmp[i])
    {
      error_log("2 dim double malloc error[2]");
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
    error_log("1 dim int malloc error");
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
    error_log("2 dim double malloc error[1]");
    exit(-1);
  }
  for(i=0; i<num_y; i++)
  {
    tmp[i] = (int *)malloc(sizeof(int)*num_x);
    if(!tmp[i])
    {
      error_log("2 dim double malloc error[2]");
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
    warning_log(msg);
    return -1;
  }else{
    return 0;
  }
}

void init_ver(char *c_matrix, 
    char *c_outer_solver, int *i_outer_maxloop, double *d_outer_eps, int *i_outer_restart, int *i_outer_kskip, int *i_outer_fix,
    char *c_inner_solver, int *i_inner_maxloop, double *d_inner_eps, int *i_inner_restart, int *i_inner_kskip, int *i_inner_fix,
   int *i_thread, bool *f_cuda, bool *f_verbose)
{
  c_matrix=NULL;
  
  c_outer_solver=NULL;
  *i_outer_maxloop=10000;
  *d_outer_eps=1e-8;
  *i_outer_restart=1000;
  *i_outer_fix=2;

  c_inner_solver=NULL;
  *i_inner_maxloop=100;
  *d_inner_eps=1e-1;
  *i_inner_restart=10;
  *i_inner_fix=2;

  *i_thread=8;
  f_cuda=false;
  f_verbose=false;
}

