/**
 * @file
 * @author Ryan Moore
 * @brief I contain the predict function for ONF.
 */

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <rya.h>

#include "onf.h"
#include "predict.h"
#include "array.h"
#include "const.h"
#include "file.h"
#include "onf.h"
#include "tommyarray.h"
#include "rlib.h"

#include "onf_math.h"

#ifdef OPENMP
#include <omp.h>
#endif

#define PSEUDO_COUNT = 16
#define TODO_VIRUS_GENOME_LENGTH = 32000


int main(int argc, char* argv[])
{
  if (argc != 5) {
    fprintf(stderr, "\n\nonf %s\n\n", ONF_VERSION);
    fprintf(stderr, "usage: %s num_threads /dir/with/host/counts /dir/with/virus/counts outdir\n",
            argv[0]);

    return 1;
  }

  char* arg_num_threads = argv[1];
  // TODO actually set the number of threads

  char* arg_host_dir   = argv[2];
  char* arg_virus_dir  = argv[3];
  char* arg_output_dir = argv[4];

  // Make the outdir if it does not exist.
  struct stat st;

  // Check for existence.
  if (stat(arg_output_dir, &st) < 0) {
    // If not, create the directory.
    // read/write/search permissions for owner and group.
    // read/search permissions for others.
    if (mkdir(arg_output_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0) {
      fprintf(stderr, "ERROR -- cannot create directory %s\n", arg_output_dir);

      return 1;
    }
  }

#ifdef OPENMP
  omp_set_num_threads(3);
#endif

  // Get count file names.

  tommy_array* host_fnames = onf_file_files_in_dir(arg_host_dir);
  assert(host_fnames);
  size_t num_input_files_host = tommy_array_size(host_fnames);

  tommy_array* virus_fnames = onf_file_files_in_dir(arg_virus_dir);
  assert(virus_fnames);
  size_t num_input_files_virus = tommy_array_size(virus_fnames);

  tommy_array* host_counts = malloc(sizeof(tommy_array));
  assert(host_counts);
  tommy_array_init(host_counts);

  tommy_array* virus_counts = malloc(sizeof(tommy_array));
  assert(virus_counts);
  tommy_array_init(virus_counts);

  for (int i = 0; i < num_input_files_host; ++i) {
    char* fname = tommy_array_get(host_fnames, i);
    assert(fname);
    // read counts
    struct onf_rya_int_array** counts = onf_read_counts2(fname);
    assert(counts != ONF_ERROR_PTR);

    struct count_array_info* count_info = malloc(sizeof(struct count_array_info));
    assert(count_info);

    count_info->fname  = fname;
    count_info->counts = counts;

    tommy_array_insert(host_counts, count_info);
  }

  // Do the same thing for the viruses.
  for (int i = 0; i < num_input_files_virus; ++i) {
    char* fname = tommy_array_get(virus_fnames, i);
    assert(fname);
    // read counts
    struct onf_rya_int_array** counts = onf_read_counts2(fname);
    assert(counts);

    struct count_array_info* count_info = malloc(sizeof(struct count_array_info));
    assert(count_info);

    count_info->fname  = fname;
    count_info->counts = counts;

    tommy_array_insert(virus_counts, count_info);
  }

  // TODO technically, I don't have to save all the viral counts.  I can just do the calculation after reading each virus file.


  int z = 0;
  assert(num_input_files_virus == tommy_array_size(virus_counts));
  assert(num_input_files_host == tommy_array_size(host_counts));

  fprintf(stdout, "vir\thost\tll\n");

#pragma omp parallel for schedule(auto) private(z)
  for (z = 0; z < num_input_files_virus; ++z) {
    struct count_array_info* virus_count_info = tommy_array_get(virus_counts, z);

    for (int i = 0; i < num_input_files_host; ++i) {
      struct count_array_info* host_count_info = tommy_array_get(host_counts, i);

      // calculate score  todo get actuall values rather than 16 and 32,000
      double log_lik = score((host_count_info->counts), (virus_count_info->counts), 16);

      fprintf(stdout, "%s\t%s\t%lf\n", (virus_count_info->fname), (host_count_info->fname), log_lik);
    }
  }
}
