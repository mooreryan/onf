#ifndef _FILE_H
#define _FILE_H

#include "tommyarray.h"
#include "cute_files.h"
#define CUTE_FILES_IMPLEMENTATION

/**
 * @brief I give the absolute path of file names in the directory recursively.
 * @param path
 * @retval Names of files in path.
 */
tommy_array* onf_file_files_in_dir(const char* path);

#endif // _FILE_H
