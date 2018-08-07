/**
 * @file
 * @author Ryan Moore
 * @brief I contain functions for the ONF library.
 */

#ifndef _ONF_H
#define _ONF_H

#include "array.h"

/**
 * @brief Encodes a nucleotide sequence as an integer array.
 *
 * @note Case insensitive.
 * @note I will encode any non-ACTG as 4.
 *
 * @param seq The nucleotide sequence to encode.
 * @param len The length of the sequence (ignoring terminating null char).
 *
 * @retval encoded_seq e.g., { 0, 0, 1, 2, 3, 4 } for "AaCTGN"
 * @retval ONF_ERROR_PTR if there are errors
 */
struct onf_int_array* onf_encode_seq(char* seq, size_t len);

/**
 * @brief Hash the encoded seq into an integer.
 *
 * Whatever the size of the encoded seq, the output will run from 0 - 2^(size-1).  For example, when size = 2:
 *
 * kmer: [0, 0], hashed: 0
 *
 * kmer: [0, 1], hashed: 1
 *
 * kmer: [0, 2], hashed: 2
 *
 * kmer: [0, 3], hashed: 3
 *
 * kmer: [1, 0], hashed: 4
 *
 * kmer: [1, 1], hashed: 5
 *
 * kmer: [1, 2], hashed: 6
 *
 * kmer: [1, 3], hashed: 7
 *
 * kmer: [2, 0], hashed: 8
 *
 * kmer: [2, 1], hashed: 9
 *
 * kmer: [2, 2], hashed: 10
 *
 * kmer: [2, 3], hashed: 11
 *
 * kmer: [3, 0], hashed: 12
 *
 * kmer: [3, 1], hashed: 13
 *
 * kmer: [3, 2], hashed: 14
 *
 * kmer: [3, 3], hashed: 15
 *
 * @param encoded_seq For example, the output of onf_hash_encoded_seq()
 * @param len Length of the encoded seq.
 *
 * @note I will return an error value if the kmer you're trying to hash doesn't have 0 - 3 for every value.  I.e., if it has any chars other than 'AaCcTtGg'.
 *
 * @retval hashed_sequence e.g., 11 for { 2, 3 }
 * @retval ONF_ERROR_INT if there are errors
 * @retval ONF_ERROR_INT if the encoded seq has non 0, 1, 2, 3 values (i.e. non A C T G )
 */
int onf_hash_int_array(struct onf_int_array* ary);

/**
 * @brief Returns an int array for counts of the correct size for counting kmers.
 * @param size Size of the kmer to count.  Must be > 0.
 *
 * @retval Array of size 2^size filled with zeros.
 * @retval ONF_ERROR_PTR if there were errors or there was a bad size.
 */
struct onf_int_array* onf_kmer_count_array_new(size_t size);

int onf_hash_lower_order_kmer(int hashed_kmer, int how_much_lower);

struct onf_int_array* onf_count_kmers(char* seq, size_t seq_len, size_t kmer_size);

#endif // _ONF_H
