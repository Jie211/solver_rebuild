#include "start.h"

int csr_start(int argc, char *argv[])
{

  int error;

  struct Parameter para;

  init_ver(&para);

  error = get_opt(argc, argv, &para);
  if(error_handle(error, "error in get_cmd")!=0)
    return -1;

  error = check_opt(&para);
  if(error_handle(error, "error in check_cmd")!=0)
    return -1;

  show_opt(&para);
  return 0;
}
