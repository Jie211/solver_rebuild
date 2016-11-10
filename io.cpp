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
  strcpy(para->list[10], "BICG");
  strcpy(para->list[11], "KSKIPBICG");

  para->isVP=false;
  para->f_cuda=false;
  para->f_verbose=false;
  strcpy(para->c_matrix, "");
  
  para->c_outer_solver=NONE;
  para->c_inner_solver=NONE;
  para->i_outer_maxloop=10000;
  para->i_inner_maxloop=100;
  para->d_outer_eps=1e-8;
  para->d_inner_eps=1e-1;
  para->i_outer_restart=1000;
  para->i_inner_restart=10;
  para->i_outer_kskip=2;
  para->i_inner_kskip=2;
  para->i_outer_fix=2;
  para->i_inner_fix=2;
  para->i_thread=8;

  strcpy(para->bx_path, "");
  strcpy(para->ptr_path, "");
  strcpy(para->col_path, "");
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
    {"Help", no_argument, NULL, 'h'},
    {0,        0,           0,      0 },
  };

  while((opt=getopt_long_only(argc, argv, "m:S:s:L:l:E:e:R:r:K:k:F:f:t:cvh", longopts, &longindex)) != -1)
  {
    switch(opt)
    {
      case 'm':
        strcpy(para->c_matrix, optarg);
        break;
      case 'S':
        error=check_solver(optarg, &para->c_outer_solver);
        if(error_handle(error, (char*)"error in outer solver name")!=0)
          return -1;
        break;
      case 's':
        error=check_solver(optarg, &para->c_inner_solver);
        if(error_handle(error, (char*)"error in inner solver name")!=0)
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
      case 'h':
        show_help();
        exit(0);
      default:
        warning_log((char*)"Unknow option");
        return -1;
    }
  }
#ifdef EBUG
  normal_log((char*)"pass get_opt");
#endif
  return 0;
}

void show_help(void)
{
  printf("Options:\n");
  printf("\t-Matrix/-m: directory name of matrix(CRS format), directory must in  \n");
  printf("\t\t|-Solver root/Solver\n");
  printf("\t\t|-Matrix/CRS/\n");
  printf("\t-OuterSolver/-S: select solver to use, if use VP method this option mean outer solver\n");
  printf("\t-InnerSolver/-s: if use VP method this option mean inner solver\n");
  printf("Option: \n"
      "\tRequire Options:\n"
      "\t\t-m/-Matrix [matrix name]-> Matrix to solve\n"
      "\t\t\t-S/-OuterSolver [solver name]-> Select Outsider solver\n"
      "\t\t\t-s/-InnerSolver [solver name]-> Select Insider solver\n"
      "\tOther Options:\n"
      "\t\t-L/-OuterLoop [int]-> maxloop for Outer solver\n"
      "\t\t-l/-InnerLoop [int]-> maxloop for Inner solver\n"
      "\t\t-E/-OuterEPS [double]-> EPS for Outer solver\n"
      "\t\t-e/-InnerEPS [double]-> EPS for Inner solver\n"
      "\t\t-R/-OuterRestart [int]-> Restart counter for Outer solver\n"
      "\t\t-r/-InnerRestart [int]-> Restart counter for Outer solver\n"
      "\t\t-K/-OuterKskip [int]-> skip num for K-skip Outer solver\n"
      "\t\t-k/-InnerKskip [int]-> skip num for K-skip Inner solver\n"
      "\t\t-F/-OuterFix [1,2]-> BugFix for K-skip Outer solver(DEBUG)\n"
      "\t\t-f/-InnerFix [1,2]-> BugFix for K-skip Inner solver(DEBUG)\n"
      "\t\t-v/-Verbose -> verbose mode\n"
      "\t\t-t/-Thread [num]-> Thread for OpenMP\n"
      "\t\t-c/-Cuda -> Cuda mode\n"
      "\t\t-h/-Help -> Help \n");
  printf("==================================================================\n");
  printf("SupportMethod: \n"
      "\tStand-along:\n"
      "\t\tcg, cr, gcr, gmres, kskipcg\n"
      "\tOuter:\n"
      "\t\tvpcg, vpcr, vpgmres\n"
      "\tInner:\n"
      "\t\tcg, cr, gcr, gmres, kskipcg\n");
  printf("==================================================================\n");

}

int check_opt(struct Parameter *para)
{
  if(strncmp(para->c_matrix, "", 128)==0)
  {
    warning_log((char*)"Must set Matrix name");
    return -1;
  }
  if(para->c_outer_solver==NONE)
  {
    warning_log((char*)"Must set a Solver name");
    return -1;
  }
  if(para->i_outer_maxloop <=0 || para->i_inner_maxloop <=0)
  {
    warning_log((char*)"Maxloop must larger than 0");
    return -1;
  }
  if(para->d_outer_eps <=0.0 || para->d_inner_eps <=0.0)
  {
    warning_log((char*)"EPS must larger than 0.0");
    return -1;
  }
  if(para->i_outer_restart <0 || para->i_inner_restart <0)
  {
    warning_log((char*)"Restart can not smaller than 0");
    return -1;
  }
  if(para->i_outer_kskip <0 || para->i_inner_kskip <0)
  {
    warning_log((char*)"Kskip can not smaller than 0");
    return -1;
  }
  if(para->c_outer_solver==VPCG || para->c_outer_solver==VPCR || para->c_outer_solver==VPGCR || para->c_outer_solver==VPGMRES)
  {
    if(para->c_inner_solver==NONE)
    {
      warning_log((char*)"Outer solver use VP method, must set a inner solver");
      return -1;
    }else if(para->c_inner_solver==VPCG || para->c_inner_solver==VPCR || para->c_inner_solver==VPGCR || para->c_inner_solver==VPGMRES)
    {
      warning_log((char*)"Outer & Inner recurrence VP method !!");
      return -1;
    }else
    {
      para->isVP=true;
    }
  }
  if(para->c_inner_solver!=NONE)
  {
    if(para->c_outer_solver!=VPCG && para->c_outer_solver!=VPCR && para->c_outer_solver!=VPGCR && para->c_outer_solver!=VPGMRES)
    {
      warning_log((char*)"Outer method was not VP, don't set inner solver");
      return -1;
    }
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
  if(para->isVP)
  {
    printf(" Inner Solver: %s\n", para->list[para->c_inner_solver]);
    printf(" Inner MaxLoop: %d\n", para->i_inner_maxloop);
    printf(" Inner EPS: %.12e\n", para->d_inner_eps);
    printf(" Inner Restart: %d\n", para->i_inner_restart);
    printf(" Inner Kskip: %d\n", para->i_inner_kskip);
    printf(" Inner Fix: %d\n", para->i_inner_fix);
    printf("****************************************\n");
  }
}

int check_solver(char *optarg, enum SolverName *solver)
{
  if(strncmp(optarg, "CG", 2)==0 || strncmp(optarg, "cg", 2)==0)
  {
    *solver=CG;
  }else if(strncmp(optarg, "CR", 2)==0 || strncmp(optarg, "cr", 2)==0)
  {
    *solver=CR;
  }else if(strncmp(optarg, "GCR", 3)==0 || strncmp(optarg, "gcr", 3)==0)
  {
    *solver=GCR;
  }else if(strncmp(optarg, "GMRES", 5)==0 || strncmp(optarg, "gmres", 5)==0)
  {
    *solver=GMRES;
  }else if(strncmp(optarg, "KSKIPCG", 7)==0 || strncmp(optarg, "kskipcg", 7)==0)
  {
    *solver=KSKIPCG;
  }else if(strncmp(optarg, "KSKIPCR", 7)==0 || strncmp(optarg, "kskipcr", 7)==0)
  {
    *solver=KSKIPCR;
  }else if(strncmp(optarg, "VPCG", 4)==0 || strncmp(optarg, "vpcg", 4)==0)
  {
    *solver=VPCG;
  }else if(strncmp(optarg, "VPCR", 4)==0 || strncmp(optarg, "vpcr", 4)==0)
  {
    *solver=VPCR;
  }else if(strncmp(optarg, "VPGCR", 5)==0 || strncmp(optarg, "vpgcr", 5)==0)
  {
    *solver=VPGCR;
  }
  else if(strncmp(optarg, "VPGMRES", 7)==0 || strncmp(optarg, "vpgmres", 7)==0)
  {
    *solver=VPGMRES;
  }
  else if(strncmp(optarg, "BICG", 4)==0 || strncmp(optarg, "bicg", 4)==0)
  {
    *solver=BICG;
  }
  else if(strncmp(optarg, "KSKIPBICG", 9)==0 || strncmp(optarg, "kskipbicg", 9)==0)
  {
    *solver=KSKIPBICG;
  }
  else{
    warning_log((char*)"not defined solver name");
    return -1;
  }
  return 0;
}

int find_mat(struct Parameter *para)
{
  DIR *dir;
  struct dirent *dp;
  struct stat st;
  char path[512]="";
  char fullpath[512+512]="";
  char search_name[512]="";

  bool dir_found=false;
  bool file_found=false;
  bool bx=false;
  bool col=false;
  bool ptr=false;

  //TODO: option to change path
  strcpy(search_name, para->c_matrix);
  strcpy(path, "../Matrix/CSR/");


  if((dir = opendir(path)) == NULL)
  {
    warning_log((char*)"open path failed");
    return -1;
  }
  
  for(dp=readdir(dir); dp!=NULL; dp=readdir(dir))
  {
    stat(dp->d_name, &st);
    if(S_ISDIR(st.st_mode))
    {
      if(strcmp(dp->d_name, search_name) == 0)
      {
        dir_found = true;
        break;
      }
    }
  }

  if(!dir_found)
  {
    warning_log((char*)"find matrix failed");
    return -1;
  }else
  {
#ifdef EBUG
    normal_log((char*)"find matrix dir");
#endif
  }

  strcpy(fullpath, path);
  strcat(fullpath, search_name);
  strcat(fullpath, "/");
  if((dir = opendir(fullpath)) == NULL)
  {
    warning_log((char*)"open fullpath failed");
    return -1;
  }
  for(dp=readdir(dir); dp!=NULL; dp=readdir(dir))
  {
    stat(dp->d_name, &st);
    /* if(S_ISDIR(st.st_mode)) */
    /* { */

      /* printf("!!%s\n", dp->d_name); */
      if(strcmp(dp->d_name, "bx.txt") == 0)
      {
        bx=true;
      }else if(strcmp(dp->d_name, "ColVal.txt") == 0)
      {
        col=true;
      }else if(strcmp(dp->d_name, "Ptr.txt") == 0)
      {
        ptr=true;
      }
      if(bx && col && ptr)
      {
        file_found=true;
        break;
      }

    /* } */
  }
  
  if(!file_found)
  {
    warning_log((char*)"matrix file meybe missing ?");
    return -1;
  }else{
    normal_log((char*)"find matrix file");
  }

  strcpy(para->bx_path, fullpath);
  strcpy(para->ptr_path, fullpath);
  strcpy(para->col_path, fullpath);
  strcat(para->bx_path, "bx.txt");
  strcat(para->ptr_path, "Ptr.txt");
  strcat(para->col_path, "ColVal.txt");

  closedir(dir);

  return 0;
}

int get_mat_head(const struct Parameter *para, int *N, int *NNZ)
{
  FILE *in1, *in2, *in3;
  int N_x[3], N_y[3], N_NNZ[3], i;

  if((in1=fopen(para->bx_path, "r")) == NULL)
  {
    warning_log((char*)"fopen failed");
    return -1;
  }
  if((in2=fopen(para->col_path, "r")) == NULL)
  {
    warning_log((char*)"fopen failed");
    return -1;
  }
  if((in3=fopen(para->ptr_path, "r")) == NULL)
  {
    warning_log((char*)"fopen failed");
    return -1;
  }

  fscanf(in1, "%d %d %d\n", &N_x[0], &N_y[0], &N_NNZ[0]);
  fscanf(in2, "%d %d %d\n", &N_x[1], &N_y[1], &N_NNZ[1]);
  fscanf(in3, "%d %d %d\n", &N_x[2], &N_y[2], &N_NNZ[2]);

  for(i=0;i<3;i++)
  {
    if(N_x[i] != N_y[i])
    {
      warning_log((char*)"N_x != N_y");
      return -1;
    }
  }
  if(N_x[0] != N_x[1] || N_x[1] != N_x[2] || N_x[2] != N_x[0])
  {
    warning_log((char*)"N_x was not same in 3files");
    return -1;
  }
  if(N_y[0] != N_y[1] || N_y[1] != N_y[2] || N_y[2] != N_y[0])
  {
    warning_log((char*)"N_y was not same in 3files");
    return -1;
  }
  if(N_NNZ[0] != N_NNZ[1] || N_NNZ[1] != N_NNZ[2] || N_NNZ[2] != N_NNZ[0])
  {
    warning_log((char*)"N_NNZ was not same in 3files");
    return -1;
  }

  *N=N_x[0];
  *NNZ=N_NNZ[0];

  fclose(in1);
  fclose(in2);
  fclose(in3);

#ifdef EBUG
  char size[10];
  sprintf(size, "get N=%d, NNZ=%d", *N, *NNZ);
  normal_log((char*)size);
#endif

  return 0;
}

int get_mat_data(const struct Parameter *para, int *col, int *ptr, double *val, double *bvec, double *xvec, const int N, const int NNZ)
{
  FILE *in1, *in2, *in3;
  int i, getint, skip1, skip2, skip3;
  double getdouble, getdouble2;

  if((in1 = fopen(para->col_path, "r")) == NULL)
  {
    warning_log((char*)"open file failed");
    return -1;
  }
  if((in2 = fopen(para->ptr_path, "r")) == NULL)
  {
    warning_log((char*)"open file failed");
    return -1;
  }
  if((in3 = fopen(para->bx_path, "r")) == NULL)
  {
    warning_log((char*)"open file failed");
    return -1;
  }

  fscanf(in1, "%d %d %d\n", &skip1, &skip2, &skip3);
  fscanf(in2, "%d %d %d\n", &skip1, &skip2, &skip3);
  fscanf(in3, "%d %d %d\n", &skip1, &skip2, &skip3);

  for(i=0;i<NNZ;i++)
  {
    fscanf(in1, "%d %le\n", &getint, &getdouble);
    col[i] = getint;
    val[i] = getdouble;
  }
  
  for(i=0;i<N+1;i++)
  {
    fscanf(in2, "%d\n", &getint);
    ptr[i] = getint;
  }

  for(i=0;i<N;i++)
  {
    fscanf(in3, "%le %le\n", &getdouble, &getdouble2);
    bvec[i] = getdouble;
    xvec[i] = getdouble2;
  }

  fclose(in1);
  fclose(in2);
  fclose(in3);
#ifdef EBUG
  normal_log((char*)"done get Mat data");
#endif

  return 0;
}
