#include <stdio.h>
#include "start.h"
#include "tools.h"

int main(int argc, char *argv[])
{
  int error = csr_start(argc, argv);
    if(error_handle(error, (char*)"error in csr_start")!=0)
    return -1;
  return 0;
}
