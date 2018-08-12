/**
 * @file
 * @author Ryan Moore
 * @brief I contain the main function for ONF.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <rya.h>

#include "main.h"
#include "array.h"
#include "const.h"
#include "file.h"
#include "onf.h"
#include "tommyarray.h"
#include "rlib.h"


#define TIMES 100000000

#ifdef OPENMP

#include <omp.h>

#endif


int main(int argc, char* argv[])
{
  if (argc != 4) {
    fprintf(stderr, "\n\nonf %s\n\n", ONF_VERSION);
    fprintf(stderr, "usage: %s num_threads /dir/with/sequence/files outdir\n", argv[0]);

    return 1;
  }

  char* arg_num_threads = argv[1];
  char* arg_input_dir   = argv[2];
  char* arg_output_dir  = argv[3];

  size_t num_input_files = 0;
  size_t z               = 0;

  // Make the modelDir if it does not exist.
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

  //#pragma omp parallel for schedule(dynamic) private(z)
//  for (z = 0; z < TIMES; ++z) { ; }

  tommy_array* fnames = onf_file_files_in_dir(arg_input_dir);
  assert(fnames);
  num_input_files = tommy_array_size(fnames);

#pragma omp parallel for schedule(auto) private(z)
  for (z = 0; z < num_input_files; ++z) {
//    printf("%2zu. %s\n", z, (char*)tommy_array_get(fnames, z));

    char* fname = tommy_array_get(fnames, z);
    if (!fname) {
      fprintf(stderr, "ERROR -- couldn't get fname for iter %zu.  Skipping.\n", z);
      continue;
    }

    // TODO ignore any non-fasta things...

    tommy_array* seqs = onf_read_seqs(fname);
    if (!seqs) {
      fprintf(stderr, "ERROR -- couldn't read seqs for '%s' (iter %zu).  Skipping.\n", fname, z);
      continue;
    }

    struct onf_rya_int_array** counts = onf_count_seq_kmers2(seqs);
    if (!counts) {
      // TODO free seqs
      fprintf(stderr, "ERROR -- counts kmers for '%s' (iter %zu).  Skipping.\n", fname, z);
      continue;
    }

    // TODO check all the rstrings
    rstring* rstring_fname = rstring_new(fname);
    rstring* indir = rfile_dirname(rstring_fname);
    rstring* extname = rfile_extname(rstring_fname);
    rstring* base = rfile_basename2(rstring_fname, extname);
    rstring* counts_fname = rstring_format("%s/%s.onf_count_dat", arg_output_dir, rstring_data(base));

    if (!counts_fname) {
      // TODO free seqs
      // TODO free counts
      fprintf(stderr, "ERROR -- couldn't format outfname '%s' (iter %zu).  Skipping.\n", fname, z);
      continue;
    }

    rya_int ret_val = onf_write_counts2(counts, rstring_data(counts_fname));
    if (ret_val != RYA_OKAY_INT) {
      // TODO free seqs
      // TODO free counts
      // TODO free counts_fname
      fprintf(stderr, "ERROR -- couldn't write counts for '%s' (iter %zu).  Skipping.\n", fname, z);
      continue;
    }

    // TODO free seqs
    // TODO free counts
    // TODO free counts_fname
  }

  return 0;
}
