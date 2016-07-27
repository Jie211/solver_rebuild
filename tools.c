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

void error_handle(int error_code, char *msg)
{
  if(error_code!=0)
  {
    error_log(msg);
  }else{
    return;
  }
}
