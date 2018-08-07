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
  char* seq = "AaNnCcTtGgXx";
  size_t len           = 12;
  int    encoded_seq[] = {0, 0, 1, 1, 2, 2, 3, 3};

  struct onf_int_array* actual = onf_encode_seq(seq, len);

  TEST_ASSERT_EQUAL(8, actual->length);
  TEST_ASSERT_EQUAL_INT_ARRAY(encoded_seq, actual->array, actual->length);

  free(actual);
}

void test___onf_encode_seq___should_NotEncodeWeirdChars(void)
{
  char* seq = "@% a_c-?";
  size_t len           = 8;
  int    encoded_seq[] = { 0, 1 };

  struct onf_int_array* actual = onf_encode_seq(seq, len);

  TEST_ASSERT_EQUAL(2, actual->length);
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

void test___onf_hash_int_array___should_HashTheSeq(void)
{
  size_t len = 2;
  struct onf_int_array* ary = onf_int_array_new(len);
  ary->array[0] = 2;
  ary->array[1] = 3;

  int expected = 11;
  int actual   = onf_hash_int_array(ary);

  TEST_ASSERT_EQUAL_INT(expected, actual);
}

void test___onf_hash_int_array___should_ReturnAnErrorCodeIfHasBadInts(void)
{
  size_t len    = 5;
  int    actual = 0;

  struct onf_int_array* ary = onf_int_array_new(len);

  ary->array[0] = 4;
  actual = onf_hash_int_array(ary);
  TEST_ASSERT_EQUAL_INT(ONF_ERROR_INT, actual);

  ary->array[0] = -1;
  actual = onf_hash_int_array(ary);
  TEST_ASSERT_EQUAL_INT(ONF_ERROR_INT, actual);
}


void test___onf_hash_int_array___should_HandleBadInput(void)
{
  struct onf_int_array* ary = NULL;
  size_t len = 2;

  int actual = onf_hash_int_array(ary);

  TEST_ASSERT_EQUAL_INT(ONF_ERROR_INT, actual);
}

void test___onf_hash_int_array___should_HashProperly(void)
{
  int size     = 3;
  int actual   = 0;
  int expected = 0;
  struct onf_int_array* ary = onf_int_array_new(size);

  for (int a = 0; a <= size; ++a) {
    ary->array[0] = a;

    for (int b = 0; b <= size; ++b) {
      ary->array[1] = b;

      for (int c = 0; c <= size; ++c) {
        ary->array[2] = c;

        actual = onf_hash_int_array(ary);

        TEST_ASSERT_EQUAL_INT(expected, actual);
        ++expected; // No side effects in a macro!
      }
    }
  }

  free(ary);
}

//////////////////////

void test___onf_kmer_count_array_new___should_ReturnNewKmerCountArray(void)
{
  size_t size        = 2;
  size_t output_size = 16;
  int    expected[]  = {
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
  };

  struct onf_int_array* actual = onf_kmer_count_array_new(size);

  TEST_ASSERT_EQUAL(output_size, actual->length);
  TEST_ASSERT_EQUAL_INT_ARRAY(expected, actual->array, output_size);

  free(actual);
}

void test___onf_kmer_count_array_new___should_ReturnErrorIfSizeIsBad(void)
{
  size_t size = 0;

  struct onf_int_array* actual = onf_kmer_count_array_new(size);

  TEST_ASSERT_EQUAL_PTR(ONF_ERROR_PTR, actual);

  free(actual);
}

//////////////////////

void test___onf_hash_lower_order_kmer___should_ReturnHashValOfLowerOrderKmer(void)
{
  int how_much_lower = 1;

  struct onf_int_array* kmer3 = onf_int_array_new(3);
  kmer3->array[0] = 0;
  kmer3->array[1] = 1;
  kmer3->array[2] = 2;

  struct onf_int_array* kmer2 = onf_int_array_new(2);
  kmer2->array[0] = 0;
  kmer2->array[1] = 1;

  int kmer3_hash = onf_hash_int_array(kmer3);

  int expected = onf_hash_int_array(kmer2);

  int actual = onf_hash_lower_order_kmer(kmer3_hash, how_much_lower);

  TEST_ASSERT_EQUAL_INT(expected, actual);
}

void test___onf_hash_lower_order_kmer___HandlesAnyLowerOrder(void)
{
  int how_much_lower = 2;

  struct onf_int_array* kmer3 = onf_int_array_new(3);
  kmer3->array[0] = 0;
  kmer3->array[1] = 1;
  kmer3->array[2] = 2;

  struct onf_int_array* kmer1 = onf_int_array_new(1);
  kmer1->array[0] = 0;

  int kmer3_hash = onf_hash_int_array(kmer3);

  int expected = onf_hash_int_array(kmer1);

  int actual = onf_hash_lower_order_kmer(kmer3_hash, how_much_lower);

  TEST_ASSERT_EQUAL_INT(expected, actual);
}

//////////////////////

void test___onf_count_kmers___should_ReturnKmerCountsForSequence(void)
{
  // Only the non-4 kmers count!!!
  // kmers: ac, ct, tg, ga, ac
  // encoded kmers are {0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 1}
  char* seq = "actNgac"; // {0, 1, 2, 3, 0, 1 }

  size_t kmer_size = 2, seq_len = 7;

  int counts[] = {
      0, 2, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1,
      1, 0, 0, 0,
  };

  struct onf_int_array* actual = onf_count_kmers(seq, seq_len, kmer_size);

  TEST_ASSERT_EQUAL(16, actual->length);
  TEST_ASSERT_EQUAL_INT_ARRAY(counts, actual->array, 16);

  free(actual);
}

void test___onf_count_kmers___should_ReturnErrorOnBadInput(void)
{
  struct onf_int_array* actual = NULL;

  // Null pointer for seq
  actual = onf_count_kmers(NULL, 10, 2);
  TEST_ASSERT_EQUAL_PTR(ONF_ERROR_PTR, actual);

  // seq_len < 1
  actual = onf_count_kmers("actg", 0, 0);
  TEST_ASSERT_EQUAL_PTR(ONF_ERROR_PTR, actual);

  // kmer_size < 1
  actual = onf_count_kmers("actg", 10, 0);
  TEST_ASSERT_EQUAL_PTR(ONF_ERROR_PTR, actual);

  // kmer_size > seq_len
  actual = onf_count_kmers("actg", 10, 20);
  TEST_ASSERT_EQUAL_PTR(ONF_ERROR_PTR, actual);
}

//////////////////////
//
//void test___onf_lower_order_counts___should_ReturnLowerOrderCountArray(void)
//{
//
//}