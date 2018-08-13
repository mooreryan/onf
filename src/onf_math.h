#ifndef _ONF_MATH_H
#define _ONF_MATH_H

#include "array.h"

/**
 * @brief These are the special big_si counts from the count2 style functions.
 * @param host_counts2
 * @param virus_counts2
 * @retval the log-likelihood of the viral counts under the model of the host counts
 */
double score(struct onf_rya_int_array** host_counts,
             struct onf_rya_int_array** virus_counts,
             double pseudo_count);

#endif // _ONF_MATH_H
