#include "const.h"
#include "file.h"
#include "cute_files.h"
#include "tommyarray.h"

#include <assert.h>

void add_fname_to_ary(cf_file_t* file, void* ary)
{
  assert(ary);
  if (file->is_reg) {
    char* fname = strdup(file->name);
    if (fname != NULL) {
      tommy_array_insert((tommy_array*) ary, fname);
    }
  }
}

tommy_array* onf_file_files_in_dir(const char* path)
{
  tommy_array* ary = malloc(sizeof(tommy_array));
  if (ary == NULL) { return ONF_ERROR_PTR; }

  tommy_array_init(ary);

  cf_traverse(path, add_fname_to_ary, (void*) ary);

  return ary;
}
