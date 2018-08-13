/**
 * @file
 * @author Ryan Moore
 * @brief I contain the build function for ONF.
 */

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <rya_file.h>
#include <rya_format.h>

#include "build.h"
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
  if (argc != 5) {
    fprintf(stderr, "\n\nonf %s\n\n", ONF_VERSION);
    fprintf(stderr, "usage: %s num_threads /dir/with/host/files /dir/with/virus/files outdir\n",
            argv[0]);

    return 1;
  }


  char* arg_num_threads = argv[1];
  // TODO actually set the number of threads

  char* arg_host_dir   = argv[2];
  char* arg_virus_dir  = argv[3];
  char* arg_output_dir = argv[4];

  size_t num_input_files       = 0;
  size_t num_input_files_virus = 0;
  size_t num_input_files_host  = 0;
  size_t z                     = 0;

  char* host_count_dir = rya_format("%s/%s", arg_output_dir, "host");
  assert(host_count_dir != RYA_ERROR_PTR);

  char* virus_count_dir = rya_format("%s/%s", arg_output_dir, "virus");
  assert(virus_count_dir != RYA_ERROR_PTR);


  // Check for existence.
  if (rya_file_exist(arg_output_dir) == rya_false) { // todo also check against rya_error
    if (mkdir(arg_output_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0) {
      // read/write/search permissions for owner and group.
      // read/search permissions for others.
      fprintf(stderr, "ERROR -- cannot create main output directory %s\n", arg_output_dir);

      return 1;
    }
  }

  // Check for host counts output dir
  if (rya_file_exist(host_count_dir) == rya_false) {
    if (mkdir(host_count_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0) {
      fprintf(stderr, "ERROR -- cannot create directory %s\n", host_count_dir);

      return 1;
    }
  }

  // Check for virus counts dir.
  if (rya_file_exist(virus_count_dir) == rya_false) {
    if (mkdir(virus_count_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0) {
      fprintf(stderr, "ERROR -- cannot create directory %s\n", virus_count_dir);

      return 1;
    }
  }


#ifdef OPENMP
  omp_set_num_threads(3);
#endif


  tommy_array* fnames = malloc(sizeof(tommy_array));
  assert(fnames);
  tommy_array_init(fnames);


  tommy_array* host_fnames = onf_file_files_in_dir(arg_host_dir);
  assert(host_fnames);
  num_input_files_host = tommy_array_size(host_fnames);

  tommy_array* virus_fnames = onf_file_files_in_dir(arg_virus_dir);
  assert(virus_fnames);
  num_input_files_virus = tommy_array_size(virus_fnames);

  num_input_files = num_input_files_host + num_input_files_virus;

  for (int i = 0; i < num_input_files_host; ++i) {
    tommy_array_insert(fnames, tommy_array_get(host_fnames, i));
  }

  for (int i = 0; i < num_input_files_virus; ++i) {
    tommy_array_insert(fnames, tommy_array_get(virus_fnames, i));
  }

  fprintf(stderr, "LOG -- counting host kmers\n");
  // Now do the hosts.
#pragma omp parallel for schedule(auto) private(z)
  for (z = 0; z < num_input_files_host; ++z) {

    char* fname = tommy_array_get(host_fnames, z);
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
    rstring* indir         = rfile_dirname(rstring_fname);
    rstring* extname       = rfile_extname(rstring_fname);
    rstring* base          = rfile_basename2(rstring_fname, extname);
    rstring* counts_fname  = rstring_format("%s/%s.onf", host_count_dir, rstring_data(base));

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

  fprintf(stderr, "LOG -- counting viral kmers\n");

  // Now do the viruses
#pragma omp parallel for schedule(auto) private(z)
  for (z = 0; z < num_input_files_virus; ++z) {

    char* fname = tommy_array_get(virus_fnames, z);
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
    rstring* indir         = rfile_dirname(rstring_fname);
    rstring* extname       = rfile_extname(rstring_fname);
    rstring* base          = rfile_basename2(rstring_fname, extname);
    rstring* counts_fname  = rstring_format("%s/%s.onf", virus_count_dir, rstring_data(base));

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
