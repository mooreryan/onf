#include "unity.h"
#include "const.h"
#include "onf.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

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
  int    encoded_seq[] = {0, 1};

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

void test___onf_hash_int_array___should_HashBigThingsFine(void)
{
  // These values are calculated with sandbox/hash_nuc.rb

  // actNgaNctNgaNctg
  int length = 12;
  int ints[] = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
  struct onf_int_array* ary = onf_int_array_new(12);
  ary->array = (int*) &ints;

  int hashed = onf_hash_int_array(ary);
  TEST_ASSERT_EQUAL(1776411, hashed);
}

void test___onf_hash_int_array___should_BigHashTest1(void)
{
  // I use signed ints, so I can handle kmers up to k length 15 without overflowing.
// ggggggggggggggg
  int length = 15;
  int ints[] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
  struct onf_int_array* ary = onf_int_array_new(15);
  ary->array = (int*) &ints;

  int hashed = onf_hash_int_array(ary);
  TEST_ASSERT_EQUAL(1073741823, hashed);
}

void test___onf_hash_int_array___should_ReturnErrorOnOverflow(void)
{
  int len    = 16; // should be big enough to overflow... 4 ** 16
  int ints[] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};

  struct onf_int_array* ary = onf_int_array_new(len);
  ary->array = (int*) &ints;

  int hashed = onf_hash_int_array(ary);

  TEST_ASSERT_EQUAL(ONF_ERROR_INT, hashed);
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

void test___onf_count_kmers___should_CountKmersUpTo15(void)
{
  // We should be able to count kmers up to 15 in length, but doing the assertion on that large of an array (like 4gb worth) is slow.
  size_t ksize = 9; // 9 is the biggest kmer we hash for big_simon app
  char* seq = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaC"; // it should be 30 long
  size_t seq_len = strlen(seq);
  assert(seq_len == 30);

  int num_kmers_in_seq     = seq_len - ksize + 1;
  int total_possible_kmers = pow(4, ksize);

  struct onf_int_array* actual_counts = onf_count_kmers(seq, seq_len, ksize);

  int* expected_counts = calloc(total_possible_kmers, sizeof(int));
  // they are all the same kmer "a" * ksize except the last one
  expected_counts[0] = num_kmers_in_seq - 1;
  // this one is ("a" * ksize-1) + "c", which will hash to 1.
  expected_counts[1] = 1;

  TEST_ASSERT_EQUAL(total_possible_kmers, actual_counts->length);
  TEST_ASSERT_EQUAL_INT_ARRAY(expected_counts, actual_counts->array, total_possible_kmers);

  onf_int_array_free(actual_counts);
  free(expected_counts);
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

//void print_int_array(struct onf_int_array* ary)
//{
//  printf("ary (%zu):", ary->length);
//  for (size_t i = 0; i < ary->length; ++i) {
////    printf(" %d", ary->array[i]);
//  }
//  putchar('\n');
//}
//
//void print_non_zero(struct onf_int_array* ary, char* msg)
//{
//  printf("%s:", msg);
//  for (size_t i = 0; i < ary->length; ++i) {
//    if (ary->array[i] != 0) {
//      printf(" %zu: %d,", i, ary->array[i]);
//    }
//  }
//  printf("\n");
//}

void test___onf_count_kmers2___should_ReturnLowerOrderCountArray(void)
{
//char* seq = "tccccccccg";
  char* seq = "actgACTGactgACTG";

  size_t seq_len = strlen(seq);

  struct onf_int_array* counts6 = onf_count_kmers(seq, seq_len, 6);
  struct onf_int_array* counts8 = onf_count_kmers(seq, seq_len, 8);
  struct onf_int_array* counts9 = onf_count_kmers(seq, seq_len, 9);
  assert(counts6->length == pow(4, 6));
  assert(counts8->length == pow(4, 8));
  assert(counts9->length == pow(4, 9));

  onf_count_kmers2(seq, seq_len);

  // The actual function being tested.
  struct onf_int_array** arrays = onf_count_kmers2(seq, seq_len);

  struct onf_int_array* actual_counts9 = arrays[2];
  struct onf_int_array* actual_counts8 = arrays[1];
  struct onf_int_array* actual_counts6 = arrays[0];

  TEST_ASSERT_EQUAL(counts9->length, actual_counts9->length);
  TEST_ASSERT_EQUAL(counts8->length, actual_counts8->length);
  TEST_ASSERT_EQUAL(counts6->length, actual_counts6->length);

  TEST_ASSERT_EQUAL_INT_ARRAY(counts9->array, actual_counts9->array, counts9->length);
  TEST_ASSERT_EQUAL_INT_ARRAY(counts8->array, actual_counts8->array, counts8->length);
  TEST_ASSERT_EQUAL_INT_ARRAY(counts6->array, actual_counts6->array, counts6->length);

  onf_int_array_free(counts6);
  onf_int_array_free(counts8);
  onf_int_array_free(counts9);

  for (int i = 0; i < 2; ++i) {
    onf_int_array_free(arrays[i]);
  }
  free(arrays);
}

void test___onf_count_kmers2___should_ReturnErrorIfSeqLenTooSmall(void)
{
  char* seq = "actg";
  size_t seq_len = 4;

  TEST_ASSERT_NULL(onf_count_kmers2(seq, seq_len));
}

void test___onf_count_kmers2___should_ReturnErrorIfSeqIsBad(void)
{
  char* seq = NULL;
  size_t seq_len = 4;

  TEST_ASSERT_NULL(onf_count_kmers2(seq, seq_len));
}