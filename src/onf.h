/**
 * @file
 * @author Ryan Moore
 * @brief I contain functions for the ONF library.
 */

#ifndef _ONF_H
#define _ONF_H

#include <stdint.h>

#include "array.h"

/**
 * @brief Encodes a nucleotide sequence as an integer array.
 *
 * @note Case insensitive.
 * @note I will skip any non-AaCcTtGg char.
 *
 * @param seq The nucleotide sequence to encode.
 * @param len The length of the sequence (ignoring terminating null char).
 *
 * @retval encoded_seq e.g., { 0, 0, 1, 2, 3 } for "AaCTGN"
 * @retval ONF_ERROR_PTR if there are errors
 */
struct onf_rya_int_array* onf_encode_seq(char* seq, size_t len);

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
 *
 * @note I will return an error value if the kmer you're trying to hash doesn't have 0 - 3 for every value.  I.e., if it has any chars other than 'AaCcTtGg'.
 *
 * @warning If you try and hash an rya_int_array with length > 15, I will return an error as you may overflow the integer!
 *
 * @retval hashed_sequence e.g., 11 for { 2, 3 }
 * @retval ONF_ERROR_INT if there are errors
 * @retval ONF_ERROR_INT if the encoded seq has non 0, 1, 2, 3 values (i.e. non A C T G ).
 * @retval ONF_ERROR_INT if rya_int_array size is > 15.  Anything bigger than this could overflow the integer.
 */
rya_int onf_hash_rya_int_array(struct onf_rya_int_array* ary);

/**
 * @brief Returns an rya_int array for counts of the correct size for counting kmers.
 * @param size Size of the kmer to count.  Must be > 0.
 *
 * @retval Array of size 2^size filled with zeros.
 * @retval ONF_ERROR_PTR if there were errors or there was a bad size.
 */
struct onf_rya_int_array* onf_kmer_count_array_new(size_t size);

/**
 * @brief Hash the lower order kmer given the hash value of the current order.
 *
 * If the kmer you hashed was { 0, 1, 2, 3 } and you set how_much_lower to 2, then I will hash { 0, 1 } without going through the whole hash calculation function.
 *
 * @note Doing it this way should be faster than calculating the whole hash for the head of the kmer seperately.
 *
 * @param hashed_kmer The hash value of a kmer.
 * @param how_much_lower How many orders lower?
 * @retval hashed kmer of a lower order
 */
rya_int onf_hash_lower_order_kmer(rya_int hashed_kmer, rya_int how_much_lower);

/**
 * @brief Count the kmers in the given seq.
 *
 * @note that the count array I return is indexed by the hash val of the kmer.
 *
 * @param seq E.g., "actg"
 * @param seq_len E.g., 4
 * @param kmer_size Size of the kmer to count
 * @retval an integer array with counts for all possible kmers of kmer_size size
 */
struct onf_rya_int_array* onf_count_kmers(char* seq, size_t seq_len, size_t kmer_size);

/**
 * @brief Like onf_count_kmers() but this one counts kmer_size 9, 8, and 6 at the same time.
 *
 * This is for the big_simon program.
 *
 * @note You should call this on all hosts.  For viruses, you only need to call onf_count_kmers() with kmer_size of 6.
 *
 * @param seq
 * @param seq_len
 * @retval An array of struct onf_rya_int_array pointers containing { counts6, counts8, counts9 }
 */
struct onf_rya_int_array** onf_count_kmers2(char* seq, size_t seq_len);

#endif // _ONF_H
