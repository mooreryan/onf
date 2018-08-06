#include "unity.h"
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
  size_t len = 12;
  int encoded_seq[] = { 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 4 };

  int* actual = onf_encode_seq(seq, len);

  TEST_ASSERT_EQUAL_INT_ARRAY(encoded_seq, actual, len);
}

void test___onf_encode_seq___should_EncodeWeirdCharsAs4(void)
{
  char* seq = "@% _-?";
  size_t len = 6;
  int encoded_seq[] = { 4, 4, 4, 4, 4, 4 };

  int* actual = onf_encode_seq(seq, len);

  TEST_ASSERT_EQUAL_INT_ARRAY(encoded_seq, actual, len);
}

void test___onf_encode_seq___should_ReturnNullWithBadString(void)
{
  char* seq = NULL;
  size_t len = 6;

  int* actual = onf_encode_seq(seq, len);

  TEST_ASSERT_EQUAL_PTR(ONF_ERROR_PTR, actual);
}

////////////////////

void test___onf_hash_encoded_seq___should_HashTheSeq(void)
{
  int encoded_seq[] = { 2, 3 };
  size_t len = 2;
  int expected = 11;
  int actual = onf_hash_encoded_seq(encoded_seq, len);

  TEST_ASSERT_EQUAL_INT(expected, actual);
}

void test___onf_hash_encoded_seq___should_HandleBadInput(void)
{
  int* encoded_seq = NULL;
  size_t len = 2;

  int actual = onf_hash_encoded_seq(encoded_seq, len);

  TEST_ASSERT_EQUAL_INT(ONF_ERROR_INT, actual);
}

void test___onf_hash_encoded_seq___should_HashProperly(void)
{
  int size = 3;
  int actual = 0;
  int expected = 0;
  int* encoded_seq = malloc(sizeof(int) * size);

  for (int a = 0; a <= size; ++a) {
    encoded_seq[0] = a;

    for (int b = 0; b <= size; ++b) {
      encoded_seq[1] = b;

      for (int c = 0; c <= size; ++c) {
        encoded_seq[2] = c;

        actual = onf_hash_encoded_seq(encoded_seq, size);

        TEST_ASSERT_EQUAL_INT(expected, actual);
        ++expected; // No side effects in a macro!
      }
    }
  }
}

////////////////////

void test___onf_kmer_count_array_new___should_ReturnNewKmerCountArray(void)
{
  size_t size = 2;
  size_t output_size = 16;
  int expected[] = {
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
  };

  int* actual = onf_kmer_count_array_new(size);

  TEST_ASSERT_EQUAL_INT_ARRAY(expected, actual, output_size);
}

void test___onf_kmer_count_array_new___should_ReturnErrorIfSizeIsBad(void)
{
  size_t size = 0;

  int* actual = onf_kmer_count_array_new(size);

  TEST_ASSERT_EQUAL_PTR(ONF_ERROR_PTR, actual);
}