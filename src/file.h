#ifndef _FILE_H
#define _FILE_H

#include <rya.h>
#include "tommyarray.h"
#include "cute_files.h"
#define CUTE_FILES_IMPLEMENTATION

/**
 * @brief I give the absolute path of file names in the directory recursively.
 * @param path
 * @retval Names of files in path.
 */
tommy_array* onf_file_files_in_dir(const char* path);

struct onf_rya_int_array** onf_read_counts2(const char* fname);

/**
 * @brief Write seq lengths.
 * @param seqs These are seq_rec structs.
 * @warning This relies on the input and output files being in the same order.
 * @return 
 */
rya_int write_seq_lengths(tommy_array* seqs);


#endif // _FILE_H
