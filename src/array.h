#ifndef _ARRAY_H
#define _ARRAY_H

#include <stdlib.h>

#include "const.h"
#include "rya.h"

#define ONF_ARRAY_DECLARATIONS(type)                                    \
  struct onf_##type##_array                                             \
  {                                                                     \
    size_t length;                                                      \
    type * array;                                                       \
  };                                                                    \
                                                                        \
  struct onf_##type##_array* onf_##type##_array_new(size_t length);     \
  void onf_##type##_array_free(struct onf_##type##_array* ary);         \
  int onf_##type##_array_bad(struct onf_##type##_array* ary);           \

#define ONF_ARRAY_FUNCTIONS(type)                                       \
  struct onf_##type##_array* onf_##type##_array_new(size_t length)      \
  {                                                                     \
    if (length < 1) { return ONF_ERROR_PTR; }                           \
                                                                        \
    struct onf_##type##_array* ary =                                    \
      malloc(sizeof(struct onf_##type##_array));                        \
                                                                        \
    if (ary == NULL) { return ONF_ERROR_PTR; }                          \
                                                                        \
    ary->length = length;                                               \
    ary->array = calloc(length, sizeof(type));                          \
                                                                        \
    /* Check mem, if NULL, need to free struct then bail. */            \
    if (ary->array == NULL) {                                           \
      free(ary);                                                        \
      return ONF_ERROR_PTR;                                             \
    }                                                                   \
                                                                        \
    return ary;                                                         \
  }                                                                     \
                                                                        \
  void onf_##type##_array_free(struct onf_##type##_array* ary)          \
  {                                                                     \
    if (ary != NULL) {                                                  \
      free(ary->array);                                                 \
    }                                                                   \
                                                                        \
    free(ary);                                                          \
  }                                                                     \
                                                                        \
  int onf_##type##_array_bad(struct onf_##type##_array* ary)            \
  {                                                                     \
    if (ary == NULL || ary->array == NULL || ary->length < 1) {         \
      return 1;                                                         \
    } else {                                                            \
      return 0;                                                         \
    }                                                                   \
  }                                                                     \

ONF_ARRAY_DECLARATIONS(rya_int)

#endif // _ARRAY_H
