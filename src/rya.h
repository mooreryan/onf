#ifndef _RYA_H
#define _RYA_H

#include <stdint.h>

#define RYA_ERROR_INT -10
#define RYA_ERROR_PTR NULL

enum rya_bool_enum { rya_false = 0, rya_true = 1 };

typedef enum rya_bool_enum rya_bool;
typedef int32_t rya_int;

rya_bool rya_file_exist(const char* fname);

#endif // _RYA_H
