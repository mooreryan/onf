#include "const.h"
#include "file.h"
#include "cute_files.h"
#include "tommyarray.h"
#include "array.h"
#include "onf.h"
#include <rya_file.h>
#include <math.h>
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


/**
 * @brief Read the big_simon counts from disc.
 * @remark We use a goto chain to clean up resources as suggested by MEM12-C of SEI CERT C Coding Standard.
 * @param fname
 * @retval ONF_ERROR_PTR on error
 * @retval counts on success
 */
struct onf_rya_int_array** onf_read_counts2(const char* fname)
{
  void* ret_val = ONF_ERROR_PTR;

  if (fname == NULL) {
    fprintf(stderr, "ERROR -- fname was null\n");

    return ONF_ERROR_PTR;
  }

  if (rya_file_exist(fname) != rya_true) {
    fprintf(stderr, "ERROR -- file doesn't exist: %s\n", fname);

    return ONF_ERROR_PTR;
  }

  FILE* infile = fopen(fname, "rb");
  if (infile == NULL) {
    fprintf(stderr, "ERROR -- couldn't open %s for reading\n", fname);

    return ONF_ERROR_PTR;
  }

  // Make the struct to hold the output.
  struct onf_rya_int_array** count_arrays = malloc(3 * sizeof(struct onf_rya_int_array*));
  if (count_arrays == NULL) {
    fprintf(stderr, "ERROR -- couldn't alloc count_arrays\n");

    goto fail_alloc_count_arrays;
  }


  int    i      = 0;
  size_t nitems = 0;

  for (i = 0; i < 3; ++i) {
    // First read number of items.
    int ary_len = 0;

    nitems = fread(&ary_len, sizeof(int), 1, infile);
    if (nitems != 1) {
      fprintf(stderr, "ERROR -- didn't read correct number of items while reading the count.  expected 1, got %zu\n",
              nitems);

      goto fail_fread_len;
    }

    // ary_len = 4 ** ksize
    int ksize = round(log(ary_len) / log(4));

    count_arrays[i] = onf_kmer_count_array_new(ksize);
    if (count_arrays[i] == NULL) {
      fprintf(stderr, "ERROR -- couldn't allocate count array iter %d\n", i);

      goto fail_alloc_array_new;

    }

    nitems = fread(count_arrays[i]->array, sizeof(*count_arrays[i]->array), count_arrays[i]->length, infile);

    if (nitems != count_arrays[i]->length) {
      // TODO do I need to check for !feof(infile) as well?
      fprintf(stderr,
              "ERROR -- didn't read correct number of items.  expected 1, got %zu\n",
              nitems);

      goto fail_fread_data;
    }
  }

  // Clean up in the reverse order of allocating.

  // If we are here, there were no errors.
  ret_val = count_arrays;
  fclose(infile);
  goto end; // skip all the clean up.

 fail_fread_data:
  // Fall to the next one.  This is just for clarity.
  ;
 fail_alloc_array_new:
  ++i; // Need to also free the i'th count_array.

  // Need to free up to the i'th count array but not the i'th one.
 fail_fread_len:
  for (int z = 0; z < i + 1; ++z) {
    onf_rya_int_array_free(count_arrays[z]);
  }
  // And to free the count_arrays.
  rya_free(count_arrays);

 fail_alloc_count_arrays:
  fclose(infile);
  ret_val = ONF_ERROR_PTR;

 end:
  return ret_val;
}

//// We treat each record in a single file as part of a single genome
//rya_int write_seq_lengths(tommy_array* seqs)
//{
//  seq_rec* rec = NULL;
//  
//  struct onf_rya_int_array* lengths = onf_rya_int_array_new(tommy_array_size(seqs));
//  
//  for (int i = 0; i < tommy_array_size(seqs); ++i) {
//    rec = tommy_array_get(seqs, i);
//    
//    lengths->array[i] = rec->seq_length;
//  }
//}