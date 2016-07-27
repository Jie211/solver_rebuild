#include <stdio.h>
#include "start.h"
#include "tools.h"

int main(int argc, char *argv[])
{
  int error = csr_start(argc, argv);
  error_handle(error, "csr_start error");
  return 0;
}
