#include <stdlib.h>
#include <math.h>

#include "onf.h"

/* See onf.h for documentation. */

int* onf_encode_seq(char* seq, size_t len)
{
  if (seq == NULL) { return ONF_ERROR_PTR; }

  int* encoded_seq = malloc(sizeof(int) * len);
  if (encoded_seq == NULL) { return ONF_ERROR_PTR; }

  for (size_t i = 0; i < len; ++i) {
    switch (seq[i]) {
      case 'A':
      case 'a':
        encoded_seq[i] = 0;
        break;

      case 'C':
      case 'c':
        encoded_seq[i] = 1;
        break;

      case 'T':
      case 't':
        encoded_seq[i] = 2;
        break;

      case 'G':
      case 'g':
        encoded_seq[i] = 3;
        break;

      default: /* Any wonky char gets converted to 4 */
        encoded_seq[i] = 4;
        break;
    }
  }

  return encoded_seq;
}

int onf_hash_encoded_seq(int* encoded_seq, size_t len)
{
  if (encoded_seq == NULL) { return ONF_ERROR_INT; }

  int hashed_val = 0;

  for (size_t i = 0; i < len; ++i) {
    hashed_val += ((1 << (2 * (len - 1 - i))) * encoded_seq[i]);
  }

  return hashed_val;
}

int* onf_kmer_count_array_new(size_t size)
{
  if (size < 1) { return ONF_ERROR_PTR; }

  int* counts = calloc(sizeof(int), pow(4, size));

  if (counts == NULL) { return ONF_ERROR_PTR; }

  return counts;
}
