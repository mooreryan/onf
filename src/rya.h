/**
 * @file
 * @author Ryan Moore
 * @brief I contain functions for the Rya library.
 */

#ifndef _RYA_H
#define _RYA_H

#include <stdint.h>

#define RYA_OKAY_INT 1
#define RYA_ERROR_INT -10
#define RYA_ERROR_PTR NULL

/**
 * @brief Return type for rya functions that seem like they should be true or false.
 *
 * Since errors are mapped to RYA_ERROR_INT you shouldn't just say `if (val)`, you need to say, `if (val == rya_true)`.  Or, you need to first check for errors, then treat it like a regular true/false value.
 */
enum rya_bool_enum { rya_false = 0, rya_true = 1, rya_error = RYA_ERROR_INT };

typedef enum rya_bool_enum rya_bool;
typedef int32_t rya_int;

/**
 * @brief I tell whether a file exists or not.
 *
 * Basically I'm a wrapper to stat().
 *
 * @param fname
 * @retval rya_error if fname is NULL
 * @retval rya_false if the file does NOT exist
 * @retval rya_true if the file exists
 */
rya_bool rya_file_exist(const char* fname);

#endif // _RYA_H
