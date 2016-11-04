#ifndef SHARE_H_INCLUDED__
#define SHARE_H_INCLUDED__

#define EBUG

#include <stdbool.h>

enum SolverName{
  CG,
  CR,
  GCR,
  GMRES,
  KSKIPCG,
  KSKIPCR,
  VPCG,
  VPCR,
  VPGCR,
  VPGMRES,
  BICG,
  NONE
};

struct Parameter{
  char list[20][128];
  bool isVP;
  bool f_cuda;
  bool f_verbose;
  char c_matrix[128];
  enum SolverName c_outer_solver;
  enum SolverName c_inner_solver;
  int i_outer_maxloop;
  int i_inner_maxloop;
  double d_outer_eps;
  double d_inner_eps;
  int i_outer_restart;
  int i_inner_restart;
  int i_outer_kskip;
  int i_inner_kskip;
  int i_outer_fix;
  int i_inner_fix;
  int i_thread;
  char bx_path[128];
  char ptr_path[128];
  char col_path[128];
};

#endif //SHARE_H_INCLUDED__

