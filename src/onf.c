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
  int val        = 0;

  for (size_t i = 0; i < len; ++i) {
    val = encoded_seq[i];

    // Only values between 0 and 3 are allowed.
    if (val < 0 || val > 3) { return ONF_ERROR_INT; }

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

int onf_hash_lower_order_kmer(int hashed_kmer, int how_much_lower)
{
  // hashed_kmer / pow(2, 2 * how_much_lower)
  return hashed_kmer >> (2 * how_much_lower);
}

int* onf_count_kmers(char* seq, size_t seq_len, size_t kmer_size)
{
  // Check args
  if (seq == NULL ||
      seq_len < 1 ||
      kmer_size < 1 ||
      kmer_size > seq_len) {

    return ONF_ERROR_PTR;
  }

  int* counts = onf_kmer_count_array_new(kmer_size);
  int hashed_kmer = 0;

  int* encoded_seq = onf_encode_seq(seq, seq_len);

  for (size_t i = 0; i < seq_len - kmer_size + 1; ++i) {
    hashed_kmer = onf_hash_encoded_seq(&encoded_seq[i], kmer_size);

    if (hashed_kmer != ONF_ERROR_INT) {
      ++counts[hashed_kmer];
    }
  }

  return counts;
}
