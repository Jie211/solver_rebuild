#include "io.h"

int get_cmd(int argc, char *argv[], 
    char *c_matrix, 
    char *c_outer_solver, int *i_outer_maxloop, double *d_outer_eps, int *i_outer_restart, int *i_outer_kskip, int *i_outer_fix,
    char *c_inner_solver, int *i_inner_maxloop, double *d_inner_eps, int *i_inner_restart, int *i_inner_kskip, int *i_inner_fix,
   int *i_thread, bool *f_cuda, bool *f_verbose)
{
  int opt;
  int longindex;
  struct option longopts[] = 
  {
    {"Matrix", required_argument, NULL, 'M'},
    {"OuterSolver", required_argument, NULL, 'S'},
    {"InnerSolver", optional_argument, NULL, 's'},
    {"OuterMaxLoop", optional_argument, NULL, 'L'},
    {"InnerMaxLoop", optional_argument, NULL, 'l'},
    {"OuterEPS", optional_argument, NULL, 'E'},
    {"InnerEPS", optional_argument, NULL, 'e'},
    {"OuterRestart", optional_argument, NULL, 'R'},
    {"InnerRestart", optional_argument, NULL, 'r'},
    {"OuterKskip", optional_argument, NULL, 'K'},
    {"InnerKskip", optional_argument, NULL, 'k'},
    {"OuterFix", optional_argument, NULL, 'F'},
    {"InnerFix", optional_argument, NULL, 'f'},
    {"Thread", optional_argument, NULL, 't'},
    {"Cuda", no_argument, NULL, 'c'},
    {"Verbose", no_argument, NULL, 'v'},
    {0,        0,           0,      0 },
  };

  while((opt=getopt_long_only(argc, argv, "m:S:s:L:l:E:e:R:r:K:k:F:f:t:cv", longopts, &longindex)) != -1)
  {
    printf("%d %s\n", longindex, longopts[longindex].name);
    switch(opt)
    {
      case 'm':
        c_matrix=optarg;
        break;
      case 'S':
        c_outer_solver=optarg;
        break;
      case 's':
        c_inner_solver=optarg;
        break;
      case 'L':
        *i_outer_maxloop=atoi(optarg);
        break;
      case 'l':
        *i_inner_maxloop=atoi(optarg);
        break;
      case 'E':
        *d_outer_eps=atof(optarg);
        break;
      case 'R':
        *i_outer_restart=atoi(optarg);
         break;
      case 'r':
        *i_inner_restart=atoi(optarg);
        break;
      case 'K':
        *i_outer_kskip=atoi(optarg);
        break;
      case 'k':
        *i_inner_kskip=atoi(optarg);
        break;
      case 'F':
        *i_outer_fix=atoi(optarg);
        break;
      case 'f':
        *i_inner_fix=atoi(optarg);
        break;
      case 't':
        *i_thread=atoi(optarg);
        break;
      case 'c':
        *f_cuda=true;
        break;
      case 'v':
        *f_verbose=true;
        break;
      default:
        error_log("unknow option");
        return -1;
    }
  }
  return 0;
}

int check_cmd(char *c_matrix, 
    char *c_outer_solver, int *i_outer_maxloop, double *d_outer_eps, int *i_outer_restart, int *i_outer_kskip, int *i_outer_fix,
    char *c_inner_solver, int *i_inner_maxloop, double *d_inner_eps, int *i_inner_restart, int *i_inner_kskip, int *i_inner_fix,
   int *i_thread, bool *f_cuda, bool *f_verbose)
{
  if(c_matrix==NULL)
  {
    error_log("Must set Matrix name");
    return -1;
  }
  if(c_outer_solver==NULL)
  {
    error_log("Must set a Solver name");
    return -1;
  }
  if(*i_outer_maxloop <=0 || *i_inner_maxloop <=0)
  {
    error_log("Maxloop must larger than 0");
    return -1;
  }
  if(*d_outer_eps <=0.0 || *d_inner_eps <=0.0)
  {
    error_log("EPS must larger than 0.0");
    return -1;
  }
  if(*i_outer_restart <0 || *i_inner_restart <0)
  {
    error_log("Restart can not smaller than 0");
    return -1;
  }
  if(*i_outer_kskip <0 || *i_inner_kskip <0)
  {
    error_log("Kskip can not smaller than 0");
    return -1;
  }
  //TODO
  //i_outer_fix & i_inner_fix

  //logic check
  return 0;
}

void show_cmd(char *c_matrix, 
    char *c_outer_solver, int *i_outer_maxloop, double *d_outer_eps, int *i_outer_restart, int *i_outer_kskip, int *i_outer_fix,
    char *c_inner_solver, int *i_inner_maxloop, double *d_inner_eps, int *i_inner_restart, int *i_inner_kskip, int *i_inner_fix,
   int *i_thread, bool *f_cuda, bool *f_verbose)
{
  printf("****************************************\n");
  printf(" Matrix: %s\n", c_matrix);
  printf(" OpenMP Thread: %d\n", *i_thread);
  printf(" Verbose: %d\n", *f_verbose);
  printf(" CUDA: %d\n", *f_cuda);
  printf("****************************************\n");
  printf(" Outer Solver: %s\n", c_outer_solver);
  printf(" Outer MaxLoop: %d\n", *i_outer_maxloop);
  printf(" Outer EPS: %.12e\n", *d_outer_eps);
  printf(" Outer Restart: %d\n", *i_outer_restart);
  printf(" Outer Kskip: %d\n", *i_outer_kskip);
  printf(" Outer Fix: %d\n", *i_outer_fix);
  printf("****************************************\n");
  printf(" Inner Solver: %s\n", c_inner_solver);
  printf(" Inner MaxLoop: %d\n", *i_inner_maxloop);
  printf(" Inner EPS: %.12e\n", *d_inner_eps);
  printf(" Inner Restart: %d\n", *i_inner_restart);
  printf(" Inner Kskip: %d\n", *i_inner_kskip);
  printf(" Inner Fix: %d\n", *i_inner_fix);
  printf("****************************************\n");
}
