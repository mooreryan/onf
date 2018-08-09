#include "rya.h"

#include <stdlib.h>
#include <sys/stat.h>

rya_bool rya_file_exist(const char* fname)
{
  if (fname == NULL) { return RYA_ERROR_INT; }

  struct stat st;

  if (stat(fname, &st) < 0) {
    return rya_false;
  }
  else {
    return rya_true;
  }
}