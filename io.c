#include "io.h"

void init_ver(struct Parameter *para)
{
  strcpy(para->list[0], "CG");
  strcpy(para->list[1], "CR");
  strcpy(para->list[2], "GCR");
  strcpy(para->list[3], "GMRES");
  strcpy(para->list[4], "KSKIPCG");
  strcpy(para->list[5], "KSKIPCR");
  strcpy(para->list[6], "VPCG");
  strcpy(para->list[7], "VPCR");
  strcpy(para->list[8], "VPGCR");
  strcpy(para->list[9], "VPGMRES");

  strcpy(para->c_matrix, "");
  
  para->c_outer_solver=NONE;
  para->i_outer_maxloop=10000;
  para->d_outer_eps=1e-8;
  para->i_outer_restart=1000;
  para->i_outer_fix=2;

  para->c_inner_solver=NONE;
  para->i_inner_maxloop=100;
  para->d_inner_eps=1e-1;
  para->i_inner_restart=10;
  para->i_inner_fix=2;

  para->i_thread=8;
  para->f_cuda=false;
  para->f_verbose=false;
}

int get_opt(int argc, char *argv[], struct Parameter *para)
{
  int opt;
  int longindex;
  int error;
  struct option longopts[] = 
  {
    {"Matrix", required_argument, NULL, 'm'},
    {"OuterSolver", required_argument, NULL, 'S'},
    {"InnerSolver", required_argument, NULL, 's'},
    {"OuterMaxLoop", required_argument, NULL, 'L'},
    {"InnerMaxLoop", required_argument, NULL, 'l'},
    {"OuterEPS", required_argument, NULL, 'E'},
    {"InnerEPS", required_argument, NULL, 'e'},
    {"OuterRestart", required_argument, NULL, 'R'},
    {"InnerRestart", required_argument, NULL, 'r'},
    {"OuterKskip", required_argument, NULL, 'K'},
    {"InnerKskip", required_argument, NULL, 'k'},
    {"OuterFix", required_argument, NULL, 'F'},
    {"InnerFix", required_argument, NULL, 'f'},
    {"Thread", required_argument, NULL, 't'},
    {"Cuda", no_argument, NULL, 'c'},
    {"Verbose", no_argument, NULL, 'v'},
    {0,        0,           0,      0 },
  };

  while((opt=getopt_long_only(argc, argv, "m:S:s:L:l:E:e:R:r:K:k:F:f:t:cv", longopts, &longindex)) != -1)
  {
    /* printf("%d %s\n", longindex, longopts[longindex].name); */
    switch(opt)
    {
      case 'm':
        strcpy(para->c_matrix, optarg);
        break;
      case 'S':
        error=check_solver(optarg, &para->c_outer_solver);
        if(error_handle(error, "error in outer solver name")!=0)
          return -1;
        break;
      case 's':
        error=check_solver(optarg, &para->c_inner_solver);
        if(error_handle(error, "error in inner solver name")!=0)
          return -1;
        break;
      case 'L':
        para->i_outer_maxloop=atoi(optarg);
        break;
      case 'l':
        para->i_inner_maxloop=atoi(optarg);
        break;
      case 'E':
        para->d_outer_eps=atof(optarg);
        break;
      case 'R':
        para->i_outer_restart=atoi(optarg);
        break;
      case 'r':
        para->i_inner_restart=atoi(optarg);
        break;
      case 'K':
        para->i_outer_kskip=atoi(optarg);
        break;
      case 'k':
        para->i_inner_kskip=atoi(optarg);
        break;
      case 'F':
        para->i_outer_fix=atoi(optarg);
        break;
      case 'f':
        para->i_inner_fix=atoi(optarg);
        break;
      case 't':
        para->i_thread=atoi(optarg);
        break;
      case 'c':
        para->f_cuda=true;
        break;
      case 'v':
        para->f_verbose=true;
        break;
      default:
        warning_log("Unknow option");
        return -1;
    }
  }
  return 0;
}

int check_opt(struct Parameter *para)
{
  if(strncmp(para->c_matrix, "", 128)==0)
  {
    warning_log("Must set Matrix name");
    return -1;
  }
  if(para->c_outer_solver==NONE)
  {
    warning_log("Must set a Solver name");
    return -1;
  }
  if(para->i_outer_maxloop <=0 || para->i_inner_maxloop <=0)
  {
    warning_log("Maxloop must larger than 0");
    return -1;
  }
  if(para->d_outer_eps <=0.0 || para->d_inner_eps <=0.0)
  {
    warning_log("EPS must larger than 0.0");
    return -1;
  }
  if(para->i_outer_restart <0 || para->i_inner_restart <0)
  {
    warning_log("Restart can not smaller than 0");
    return -1;
  }
  if(para->i_outer_kskip <0 || para->i_inner_kskip <0)
  {
    warning_log("Kskip can not smaller than 0");
    return -1;
  }
  return 0;
}

void show_opt(struct Parameter *para)
{
  printf("****************************************\n");
  printf(" Matrix: %s\n", para->c_matrix);
  printf(" OpenMP Thread: %d\n", para->i_thread);
  printf(" Verbose: %d\n", para->f_verbose);
  printf(" CUDA: %d\n", para->f_cuda);
  printf("****************************************\n");
  printf(" Outer Solver: %s\n", para->list[para->c_outer_solver]);
  printf(" Outer MaxLoop: %d\n", para->i_outer_maxloop);
  printf(" Outer EPS: %.12e\n", para->d_outer_eps);
  printf(" Outer Restart: %d\n", para->i_outer_restart);
  printf(" Outer Kskip: %d\n", para->i_outer_kskip);
  printf(" Outer Fix: %d\n", para->i_outer_fix);
  printf("****************************************\n");
  printf(" Inner Solver: %s\n", para->list[para->c_inner_solver]);
  printf(" Inner MaxLoop: %d\n", para->i_inner_maxloop);
  printf(" Inner EPS: %.12e\n", para->d_inner_eps);
  printf(" Inner Restart: %d\n", para->i_inner_restart);
  printf(" Inner Kskip: %d\n", para->i_inner_kskip);
  printf(" Inner Fix: %d\n", para->i_inner_fix);
  printf("****************************************\n");
}

int check_solver(char *optarg, enum SolverName *solver)
{
  if(strncmp(optarg, "CG", 2)==0 || strncmp(optarg, "cg", 2)==0)
  {
    printf("get cg\n");
    *solver=CG;
  }else if(strncmp(optarg, "CR", 2)==0 || strncmp(optarg, "cr", 2)==0)
  {
    printf("get cr\n");
    *solver=CR;
  }else if(strncmp(optarg, "GCR", 3)==0 || strncmp(optarg, "gcr", 3)==0)
  {
    printf("get gcr\n");
    *solver=GCR;
  }else if(strncmp(optarg, "GMRES", 5)==0 || strncmp(optarg, "gmres", 5)==0)
  {
    printf("get gmres\n");
    *solver=GMRES;
  }else if(strncmp(optarg, "KSKIPCG", 7)==0 || strncmp(optarg, "kskipcg", 7)==0)
  {
    printf("get kskipcg\n");
    *solver=KSKIPCG;
  }else if(strncmp(optarg, "KSKIPCR", 7)==0 || strncmp(optarg, "kskipcr", 7)==0)
  {
    printf("get kskipcr\n");
    *solver=KSKIPCR;
  }else if(strncmp(optarg, "VPCG", 4)==0 || strncmp(optarg, "vpcg", 4)==0)
  {
    printf("get vpcg\n");
    *solver=VPCG;
  }else if(strncmp(optarg, "VPCR", 4)==0 || strncmp(optarg, "vpcr", 4)==0)
  {
    printf("get vpcr\n");
    *solver=VPCR;
  }else if(strncmp(optarg, "VPGMRES", 7)==0 || strncmp(optarg, "vpgmres", 7)==0)
  {
    printf("get vpgmres\n");
    *solver=VPGMRES;
  }else{
    warning_log("not defined solver name");
    return -1;
  }
  return 0;
}
