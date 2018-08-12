#include "const.h"
#include "file.h"
#include "cute_files.h"
#include "tommyarray.h"
#include <rya_file.h>

#include <assert.h>

/**
 * @brief Callback for cf_traverse()
 *
 * I add filenames to the tommy array.
 *
 * @warning File names are malloc'd in this function, so you'll need to free them later.
 *
 * @todo Probably want to use something other than cute files' max lengths.
 *
 * @param file
 * @param ary
 */
static void add_fname_to_ary(cf_file_t* file, void* ary)
{
  assert(ary);
  if (file->is_reg) {

    char* fname = strdup(file->path);
    if (fname != NULL) {
      tommy_array_insert((tommy_array*)ary, fname);
    }
  }
}

/**
 * @brief Puts the paths of all regular files in a dir into an array.
 * @note I will search recursively in the directory.
 * @warning The caller must free the returned array.
 * @param path The path to search.
 * @retval An array with the names of the regular files in the directory.
 * @retval ONF_ERROR_PTR if path is NULL or doesn't exist
 * @retval ONF_ERROR_PTR if errors
 */
tommy_array* onf_file_files_in_dir(const char* path)
{
  if (path == NULL || rya_file_exist(path) != rya_true) { return ONF_ERROR_PTR; }

  tommy_array* ary = malloc(sizeof(tommy_array));
  if (ary == NULL) { return ONF_ERROR_PTR; }

  tommy_array_init(ary);

  cf_traverse(path, add_fname_to_ary, (void*)ary);

  return ary;
}
