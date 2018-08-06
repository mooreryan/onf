#include "unity.h"
#include "const.h"
#include "onf.h"

#include <stdlib.h>

void setUp(void)
{
}

void tearDown(void)
{
}

////////////////////

void test___onf_encode_seq___should_EncodeTheSequence(void)
{
  char* seq = "AaCcTtGgNnXx";
  size_t len           = 12;
  int    encoded_seq[] = { 0, 0, 1, 1, 2, 2, 3, 3 };

  struct onf_int_array* actual = onf_encode_seq(seq, len);

  TEST_ASSERT_EQUAL(8, actual->length);
  TEST_ASSERT_EQUAL_INT_ARRAY(encoded_seq, actual->array, actual->length);

  free(actual);
}

void test___onf_encode_seq___should_NotEncodeWeirdChars(void)
{
  char* seq = "@% a_-?";
  size_t len           = 7;
  int    encoded_seq[] = { 0 };

  struct onf_int_array* actual = onf_encode_seq(seq, len);

  TEST_ASSERT_EQUAL(1, actual->length);
  TEST_ASSERT_EQUAL_INT_ARRAY(encoded_seq, actual->array, actual->length);

  free(actual);
}

void test___onf_encode_seq___should_ReturnNullWithBadString(void)
{
  char* seq = NULL;
  size_t len = 6;

  struct onf_int_array* actual = onf_encode_seq(seq, len);

  TEST_ASSERT_EQUAL_PTR(ONF_ERROR_PTR, actual);

  free(actual);
}

////////////////////

//void test___onf_hash_encoded_seq___should_HashTheSeq(void)
//{
//  size_t len           = 2;
//  // Will give {2, 3}
//  struct onf_int_array* encoded_seq = onf_encode_seq("TG", len);
//
//  int    expected      = 11;
//  int    actual        = onf_hash_encoded_seq(encoded_seq, len);
//
//  TEST_ASSERT_EQUAL_INT(expected, actual);
//}

//void test___onf_hash_encoded_seq___should_ReturnAnErrorCodeIfHasAmbigChars(void)
//{
//  int    encoded_seq[] = {0, 1, 2, 3, 4};
//  size_t len           = 5;
//  int    actual        = onf_hash_encoded_seq(encoded_seq, len);
//
//  TEST_ASSERT_EQUAL_INT(ONF_ERROR_INT, actual);
//}
//
//
//void test___onf_hash_encoded_seq___should_HandleBadInput(void)
//{
//  struct onf_int_array* encoded_seq = NULL;
//  size_t len = 2;
//
//  int actual = onf_hash_encoded_seq(encoded_seq, len);
//
//  TEST_ASSERT_EQUAL_INT(ONF_ERROR_INT, actual);
//}
//
//void test___onf_hash_encoded_seq___should_HashProperly(void)
//{
//  int size     = 3;
//  int actual   = 0;
//  int expected = 0;
//  struct onf_int_array* encoded_seq = malloc(sizeof(int) * size);
//
//  for (int a = 0; a <= size; ++a) {
//    encoded_seq[0] = a;
//
//    for (int b = 0; b <= size; ++b) {
//      encoded_seq[1] = b;
//
//      for (int c = 0; c <= size; ++c) {
//        encoded_seq[2] = c;
//
//        actual = onf_hash_encoded_seq(encoded_seq, size);
//
//        TEST_ASSERT_EQUAL_INT(expected, actual);
//        ++expected; // No side effects in a macro!
//      }
//    }
//  }
//
//  free(encoded_seq);
//}
//
//////////////////////
//
//void test___onf_kmer_count_array_new___should_ReturnNewKmerCountArray(void)
//{
//  size_t size        = 2;
//  size_t output_size = 16;
//  int    expected[]  = {
//      0, 0, 0, 0,
//      0, 0, 0, 0,
//      0, 0, 0, 0,
//      0, 0, 0, 0,
//  };
//
//  struct onf_int_array* actual = onf_kmer_count_array_new(size);
//
//  TEST_ASSERT_EQUAL_INT_ARRAY(expected, actual, output_size);
//
//  free(actual);
//}
//
//void test___onf_kmer_count_array_new___should_ReturnErrorIfSizeIsBad(void)
//{
//  size_t size = 0;
//
//  struct onf_int_array* actual = onf_kmer_count_array_new(size);
//
//  TEST_ASSERT_EQUAL_PTR(ONF_ERROR_PTR, actual);
//
//  free(actual);
//}
//
//////////////////////
//
//void test___onf_hash_lower_order_kmer___should_ReturnHashValOfLowerOrderKmer(void)
//{
//  int kmer3[]        = {0, 1, 2};
//  int kmer2[]        = {0, 1};
//  int how_much_lower = 1;
//
//  int kmer3_hash = onf_hash_encoded_seq(kmer3, 3);
//
//  int expected = onf_hash_encoded_seq(kmer2, 2);
//
//  int actual = onf_hash_lower_order_kmer(kmer3_hash, how_much_lower);
//
//  TEST_ASSERT_EQUAL_INT(expected, actual);
//}
//
//void test___onf_hash_lower_order_kmer___HandlesAnyLowerOrder(void)
//{
//  int kmer3[]        = {0, 1, 2};
//  int kmer1[]        = {0};
//  int how_much_lower = 2;
//
//  int kmer3_hash = onf_hash_encoded_seq(kmer3, 3);
//
//  int expected = onf_hash_encoded_seq(kmer1, 1);
//
//  int actual = onf_hash_lower_order_kmer(kmer3_hash, how_much_lower);
//
//  TEST_ASSERT_EQUAL_INT(expected, actual);
//}
//
//////////////////////
//
//void test___onf_count_kmers___should_ReturnKmerCountsForSequence(void)
//{
//  // encoded kmers are {0, 1}, {1, 2}, {2, 4}, {4, 3}, {3, 0}, {0, 1}
//  char* seq = "actNgac"; // {0, 1, 2, 4, 3, 0, 1}
//  // Only the non-4 kmers count!!!
//
//  size_t kmer_size = 2, seq_len = 7;
//
//  int counts[] = {
//      0, 2, 0, 0,
//      0, 0, 1, 0,
//      0, 0, 0, 0,
//      1, 0, 0, 0,
//  };
//
//  struct onf_int_array* actual = onf_count_kmers(seq, seq_len, kmer_size);
//
//  TEST_ASSERT_EQUAL_INT_ARRAY(counts, actual, 16);
//
//  free(actual);
//}
//
//void test___onf_count_kmers___should_ReturnErrorOnBadInput(void)
//{
//  struct onf_int_array* actual = NULL;
//
//  // Null pointer for seq
//  actual = onf_count_kmers(NULL, 10, 2);
//  TEST_ASSERT_EQUAL_PTR(ONF_ERROR_PTR, actual);
//
//  // seq_len < 1
//  actual = onf_count_kmers("actg", 0, 0);
//  TEST_ASSERT_EQUAL_PTR(ONF_ERROR_PTR, actual);
//
//  // kmer_size < 1
//  actual = onf_count_kmers("actg", 10, 0);
//  TEST_ASSERT_EQUAL_PTR(ONF_ERROR_PTR, actual);
//
//  // kmer_size > seq_len
//  actual = onf_count_kmers("actg", 10, 20);
//  TEST_ASSERT_EQUAL_PTR(ONF_ERROR_PTR, actual);
//}
//
//////////////////////
//
//void test___onf_lower_order_counts___should_ReturnLowerOrderCountArray(void)
//{
//
//}