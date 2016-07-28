#include "start.h"

int csr_start(int argc, char *argv[])
{

  int c_size = 100;
  bool f_cuda;
  bool f_verbose;
  /* char *c_matrix=NULL; */
  /* char *c_outer_solver=NULL; */
  /* char *c_inner_solver=NULL; */
  char c_matrix[c_size];
  char c_outer_solver[c_size];
  char c_inner_solver[c_size];
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

  int error;

  init_ver(c_matrix, 
      c_outer_solver, &i_outer_maxloop, &d_outer_eps, &i_outer_restart, &i_outer_kskip, &i_outer_fix,
      c_inner_solver, &i_inner_maxloop, &d_inner_eps, &i_inner_restart, &i_inner_kskip, &i_inner_fix,
      &i_thread, &f_cuda, &f_verbose);

  error = get_cmd(argc, argv,
      c_matrix, 
      c_outer_solver, &i_outer_maxloop, &d_outer_eps, &i_outer_restart, &i_outer_kskip, &i_outer_fix,
      c_inner_solver, &i_inner_maxloop, &d_inner_eps, &i_inner_restart, &i_inner_kskip, &i_inner_fix,
      &i_thread, &f_cuda, &f_verbose);

  /* if(error!=0) */
  /* { */
  /*   error_log("error in get_cmd"); */
  /*   return -1; */
  /* } */
  if(error_handle(error, "error in get_cmd")!=0)
    return -1;

  error = check_cmd(c_matrix, 
      c_outer_solver, &i_outer_maxloop, &d_outer_eps, &i_outer_restart, &i_outer_kskip, &i_outer_fix,
      c_inner_solver, &i_inner_maxloop, &d_inner_eps, &i_inner_restart, &i_inner_kskip, &i_inner_fix,
      &i_thread, &f_cuda, &f_verbose, c_size);

  /* if(error!=0) */
  /* { */
  /*   error_log("error in check_cmd"); */
  /*   return -1; */
  /* } */
  if(error_handle(error, "error in check_cmd")!=0)
    return -1;

  show_cmd(c_matrix, 
      c_outer_solver, &i_outer_maxloop, &d_outer_eps, &i_outer_restart, &i_outer_kskip, &i_outer_fix,
      c_inner_solver, &i_inner_maxloop, &d_inner_eps, &i_inner_restart, &i_inner_kskip, &i_inner_fix,
      &i_thread, &f_cuda, &f_verbose);

  return 0;
}
