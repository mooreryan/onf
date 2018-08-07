#include <stdlib.h>
#include <math.h>

#include <stdio.h> // For debug

#include "const.h"
#include "onf.h"

ONF_ARRAY_FUNCTIONS(int)

/* See onf.h for documentation. */

struct onf_int_array* onf_encode_seq(char* seq, size_t len)
{
  if (seq == NULL) { return ONF_ERROR_PTR; }

  struct onf_int_array* encoded_seq = onf_int_array_new(len);

  if (onf_int_array_bad(encoded_seq)) { return ONF_ERROR_PTR; }

  size_t bad_chars = 0;
  size_t new_i     = 0;

  for (size_t i = 0; i < len; ++i) {
    switch (seq[i]) {
      case 'A':
      case 'a':
        encoded_seq->array[new_i++] = 0;
        break;

      case 'C':
      case 'c':
        encoded_seq->array[new_i++] = 1;
        break;

      case 'T':
      case 't':
        encoded_seq->array[new_i++] = 2;
        break;

      case 'G':
      case 'g':
        encoded_seq->array[new_i++] = 3;
        break;

      default: /* Any wonky char gets dropped */
        ++bad_chars;
        break;
    }
  }

  // Subtract off the number of bad chars from the length.
  encoded_seq->length -= bad_chars;

  return encoded_seq;
}

int onf_hash_int_array(struct onf_int_array* ary)
{
  if (onf_int_array_bad(ary)) { return ONF_ERROR_INT; }

  int hashed_val = 0;
  int val        = 0;

  for (size_t i = 0; i < ary->length; ++i) {
    val = ary->array[i];

    // Only values between 0 and 3 are allowed.
    if (val < 0 || val > 3) { return ONF_ERROR_INT; }

    hashed_val += ((1 << (2 * (ary->length - 1 - i))) * ary->array[i]);
  }

  return hashed_val;
}

struct onf_int_array* onf_kmer_count_array_new(size_t size)
{
  if (size < 1) { return ONF_ERROR_PTR; }

  struct onf_int_array* counts = onf_int_array_new((size_t) pow(4, size));

  if (counts == NULL) { return ONF_ERROR_PTR; }

  return counts;
}

int onf_hash_lower_order_kmer(int hashed_kmer, int how_much_lower)
{
  // hashed_kmer / pow(2, 2 * how_much_lower)
  return hashed_kmer >> (2 * how_much_lower);
}

struct onf_int_array* onf_count_kmers(char* seq, size_t seq_len, size_t kmer_size)
{
  // Check args
  if (seq == NULL ||
      seq_len < 1 ||
      kmer_size < 1 ||
      kmer_size > seq_len) {

    return ONF_ERROR_PTR;
  }

  struct onf_int_array* counts = onf_kmer_count_array_new(kmer_size);
  int hashed_kmer = 0;

  struct onf_int_array* encoded_seq = onf_encode_seq(seq, seq_len);

  struct onf_int_array* tmp_ary = onf_int_array_new(encoded_seq->length);
  tmp_ary->length = kmer_size;

  for (size_t i = 0; i < encoded_seq->length - kmer_size + 1; ++i) {
    tmp_ary->array = &encoded_seq->array[i];

    hashed_kmer = onf_hash_int_array(tmp_ary);

    if (hashed_kmer != ONF_ERROR_INT) {
      ++(counts->array[hashed_kmer]);
    }
  }

  return counts;
}

void onf_count_kmers2(char* seq, size_t length, struct onf_int_array* counts6, struct onf_int_array* counts8,
                      struct onf_int_array* counts9)
{
  // Check args
  if (seq == NULL ||
      length < 9) {

//    return ONF_ERROR_PTR;
  }

  size_t kmer_size = 9, i = 0;

  int hashed_kmer9 = 0, hashed_kmer8 = 0, hashed_kmer6 = 0;

  struct onf_int_array* encoded_seq = onf_encode_seq(seq, length);

  struct onf_int_array* tmp_ary = onf_int_array_new(encoded_seq->length);
  tmp_ary->length = kmer_size;


  int num_9mers = encoded_seq->length - kmer_size + 1;

  for (i = 0; i < num_9mers; ++i) {


    tmp_ary->array = &encoded_seq->array[i];

    hashed_kmer9 = onf_hash_int_array(tmp_ary);

    hashed_kmer8 = onf_hash_lower_order_kmer(hashed_kmer9, 1);
    hashed_kmer6 = onf_hash_lower_order_kmer(hashed_kmer9, 3);

    if (hashed_kmer9 != ONF_ERROR_INT) {
      ++(counts9->array[hashed_kmer9]);
    }

    if (hashed_kmer8 != ONF_ERROR_INT) {
      ++(counts8->array[hashed_kmer8]);
    }

    if (hashed_kmer6 != ONF_ERROR_INT) {
      ++(counts6->array[hashed_kmer6]);
    }
  }

  // need to catch the last couple of kmers
  kmer_size = 8;
  tmp_ary->array  = &encoded_seq->array[i++]; // for the last 8mer
  tmp_ary->length = kmer_size;

  // Get the 8mer at this pos.
  hashed_kmer8 = onf_hash_int_array(tmp_ary);
  if (hashed_kmer8 != ONF_ERROR_INT) {
    ++(counts8->array[hashed_kmer8]);
  }

  // Get the 6mer at this pos.
  hashed_kmer6 = onf_hash_lower_order_kmer(hashed_kmer8, 2);
  if (hashed_kmer6 != ONF_ERROR_INT) {
    ++(counts6->array[hashed_kmer6]);
  }

  // Kmer size 6 to finish it off.
  kmer_size = 6;
  tmp_ary->length = kmer_size;

  int num_6mers = encoded_seq->length - kmer_size + 1;

  for (size_t z = i; z < num_6mers; ++z) {
    tmp_ary->array = &encoded_seq->array[z];
    hashed_kmer6 = onf_hash_int_array(tmp_ary);
    if (hashed_kmer6 != ONF_ERROR_INT) {
      ++(counts6->array[hashed_kmer6]);
    }
  }
}
