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
    /* printf("%d %s\n", longindex, longopts[longindex].name); */
    switch(opt)
    {
      case 'M':
        c_matrix=optarg;
        break;
      case 'S':
        c_outer_solver=optarg;
        break;
      case 's':
        c_inner_solver=optarg;
        break;
      case 'L':
        i_outer_maxloop=atoi(optarg);
        break;
      case 'l':
        i_inner_maxloop=atoi(optarg);
        break;
      case 'E':
        d_outer_eps=atof(optarg);

    }
  }
}
