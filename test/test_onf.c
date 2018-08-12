#include "unity.h"
#include "const.h"
#include "onf.h"
#include "tommyarray.h"
#include "rlib.h"
#include "rya.h"
#include "array.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <zlib.h> // For gzFile

#include <math.h> // For pow


void setUp(void)
{
}

void tearDown(void)
{
}

////////////////////

void test___onf_read_seqs___should_ReturnTheSeqs(void)
{
  char* this_fname = __FILE__;

  rstring* this_file = rstring_new(this_fname);
  assert(this_file);

  rstring* dirname = rfile_dirname(this_file);
  assert(dirname);

  rstring* path = rstring_format("%s/test_files/seqs/s1.fa", rstring_data(dirname));
  assert(path);

  tommy_array* seqs = onf_read_seqs(rstring_data(path));

  int num_seqs = 2;
  int seq_len  = 8;
  int id_len   = 2;

  seq_rec* rec = NULL;

  TEST_ASSERT_EQUAL(num_seqs, tommy_array_size(seqs));

  rec = tommy_array_get(seqs, 0);
  TEST_ASSERT_EQUAL(seq_len, rec->seq_length);
  TEST_ASSERT_EQUAL(id_len, rec->id_length);
  TEST_ASSERT_EQUAL_STRING("s1", rec->id);
  TEST_ASSERT_EQUAL_STRING("ACTGactg", rec->seq);
  seq_rec_free(rec);

  rec = tommy_array_get(seqs, 1);
  TEST_ASSERT_EQUAL(seq_len, rec->seq_length);
  TEST_ASSERT_EQUAL(id_len, rec->id_length);
  TEST_ASSERT_EQUAL_STRING("s2", rec->id);
  TEST_ASSERT_EQUAL_STRING("actgACTG", rec->seq);
  seq_rec_free(rec);

  tommy_array_done(seqs);
  free(seqs);
}

void test___onf_read_seqs___return_ErrorBadFname(void)
{
  TEST_ASSERT_EQUAL_PTR(ONF_ERROR_PTR, onf_read_seqs(NULL));
}

void test___onf_read_seqs___return_ErrorOnFileNotExist(void)
{
  TEST_ASSERT_EQUAL_PTR(ONF_ERROR_PTR, onf_read_seqs("arsotienarsoteinarsotn"));
}

////////////////////

void test___onf_count_seq_kmers2(void)
{
  rstring* this_file = rstring_new(__FILE__);
  assert(this_file);

  rstring* dirname = rfile_dirname(this_file);
  assert(dirname);

  rstring* path = rstring_format("%s/test_files/seqs/s2.fa", rstring_data(dirname));
  assert(path);

  tommy_array* seqs = onf_read_seqs(rstring_data(path));
  assert(seqs);

  struct onf_rya_int_array** counts = onf_count_seq_kmers2(seqs);
  assert(counts);

  struct onf_rya_int_array* counts6 = counts[0];
  assert(counts6);

  struct onf_rya_int_array* counts8 = counts[1];
  assert(counts8);

  struct onf_rya_int_array* counts9 = counts[2];
  assert(counts9);

  // The seqs are both "a" * 9.  And there are two of them.
  for (size_t i = 0; i < counts9->length; ++i) {
    if (i == 0) {
      TEST_ASSERT_EQUAL(1 * 2, counts9->array[0]);
      TEST_ASSERT_EQUAL(2 * 2, counts8->array[0]);
      TEST_ASSERT_EQUAL(4 * 2, counts6->array[0]);
    }
    else {
      TEST_ASSERT_EQUAL(0, counts9->array[i]);

      if (i < counts6->length) {
        TEST_ASSERT_EQUAL(0, counts6->array[i]);
      }

      if (i < counts8->length) {
        TEST_ASSERT_EQUAL(0, counts8->array[i]);
      }
    }
  }

  onf_rya_int_array_free(counts6);
  onf_rya_int_array_free(counts8);
  onf_rya_int_array_free(counts9);
  free(counts);

  rstring_free(this_file);
  rstring_free(dirname);
  rstring_free(path);
}

////////////////////

void test___onf_write_counts___should_WriteTheCounts(void)
{
  struct onf_rya_int_array* counts = onf_rya_int_array_new(3);
  assert(counts);

  counts->array[0] = 1;
  counts->array[1] = 2;
  counts->array[2] = 3;

  FILE* fstream = fopen("test.dat", "wb");
  assert(fstream);

  size_t num_written = onf_write_counts(counts, fstream);

  fclose(fstream);

  TEST_ASSERT_EQUAL(RYA_OKAY_INT, num_written);
}

////////////////////

void test___onf_write_counts2___should_WriteTheCounts(void)
{
  rstring* this_file = rstring_new(__FILE__);
  assert(this_file);

  rstring* dirname = rfile_dirname(this_file);
  assert(dirname);

  rstring* path = rstring_format("%s/test_files/seqs/s2.fa", rstring_data(dirname));
  assert(path);

  tommy_array* seqs = onf_read_seqs(rstring_data(path));
  assert(seqs);

  struct onf_rya_int_array** counts = onf_count_seq_kmers2(seqs);
  assert(counts);

  rya_int ret_val = onf_write_counts2(counts, "test.dat");

  TEST_ASSERT_EQUAL(RYA_OKAY_INT, ret_val);

  //TODO perhaps test the actual values.
}
////////////////////

void test___onf_read_counts___should_ReadTheCounts(void)
{
  struct onf_rya_int_array* counts = onf_rya_int_array_new(3);
  assert(counts);

  counts->array[0] = 1;
  counts->array[1] = 2;
  counts->array[2] = 3;

  FILE* fstream = fopen("test.dat", "wb");
  assert(fstream);

  size_t num_written = onf_write_counts(counts, fstream);
  assert(num_written == RYA_OKAY_INT);

  fclose(fstream);

  struct onf_rya_int_array* actual_counts = onf_read_counts("test.dat");
  assert(actual_counts);

  TEST_ASSERT_EQUAL(3, actual_counts->length);

  TEST_ASSERT_EQUAL_INT32_ARRAY(counts->array, actual_counts->array, 3);
}

////////////////////

void test___onf_encode_seq___should_EncodeTheSequence(void)
{
  char* seq = "AaNnCcTtGgXx";
  size_t  len           = 12;
  rya_int encoded_seq[] = {0, 0, 1, 1, 2, 2, 3, 3};

  struct onf_rya_int_array* actual = onf_encode_seq(seq, len);

  TEST_ASSERT_EQUAL(8, actual->length);
  TEST_ASSERT_EQUAL_INT_ARRAY(encoded_seq, actual->array, actual->length);

  free(actual);
}

void test___onf_encode_seq___should_NotEncodeWeirdChars(void)
{
  char* seq = "@% a_c-?";
  size_t  len           = 8;
  rya_int encoded_seq[] = {0, 1};

  struct onf_rya_int_array* actual = onf_encode_seq(seq, len);

  TEST_ASSERT_EQUAL(2, actual->length);
  TEST_ASSERT_EQUAL_INT_ARRAY(encoded_seq, actual->array, actual->length);

  free(actual);
}

void test___onf_encode_seq___should_ReturnNullWithBadString(void)
{
  char* seq = NULL;
  size_t len = 6;

  struct onf_rya_int_array* actual = onf_encode_seq(seq, len);

  TEST_ASSERT_EQUAL_PTR(ONF_ERROR_PTR, actual);

  free(actual);
}

////////////////////

void test___onf_hash_rya_int_array___should_HashTheSeq(void)
{
  size_t len = 2;
  struct onf_rya_int_array* ary = onf_rya_int_array_new(len);
  ary->array[0] = 2;
  ary->array[1] = 3;

  rya_int expected = 11;
  rya_int actual   = onf_hash_rya_int_array(ary);

  TEST_ASSERT_EQUAL_INT(expected, actual);
}

void test___onf_hash_rya_int_array___should_ReturnAnErrorCodeIfHasBadInts(void)
{
  size_t  len    = 5;
  rya_int actual = 0;

  struct onf_rya_int_array* ary = onf_rya_int_array_new(len);

  ary->array[0] = 4;
  actual = onf_hash_rya_int_array(ary);
  TEST_ASSERT_EQUAL_INT(ONF_ERROR_INT32_T, actual);

  ary->array[0] = -1;
  actual = onf_hash_rya_int_array(ary);
  TEST_ASSERT_EQUAL_INT(ONF_ERROR_INT32_T, actual);
}


void test___onf_hash_rya_int_array___should_HandleBadInput(void)
{
  struct onf_rya_int_array* ary = NULL;
  size_t len = 2;

  rya_int actual = onf_hash_rya_int_array(ary);

  TEST_ASSERT_EQUAL_INT(ONF_ERROR_INT32_T, actual);
}

void test___onf_hash_rya_int_array___should_HashProperly(void)
{
  rya_int size     = 3;
  rya_int actual   = 0;
  rya_int expected = 0;
  struct onf_rya_int_array* ary = onf_rya_int_array_new(size);

  for (rya_int a = 0; a <= size; ++a) {
    ary->array[0] = a;

    for (rya_int b = 0; b <= size; ++b) {
      ary->array[1] = b;

      for (rya_int c = 0; c <= size; ++c) {
        ary->array[2] = c;

        actual = onf_hash_rya_int_array(ary);

        TEST_ASSERT_EQUAL_INT(expected, actual);
        ++expected; // No side effects in a macro!
      }
    }
  }

  free(ary);
}

void test___onf_hash_rya_int_array___should_HashBigThingsFine(void)
{
  // These values are calculated with sandbox/hash_nuc.rb

  // actNgaNctNgaNctg
  rya_int length = 12;
  rya_int ints[] = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
  struct onf_rya_int_array* ary = onf_rya_int_array_new(12);
  ary->array = (rya_int * ) & ints;

  rya_int hashed = onf_hash_rya_int_array(ary);
  TEST_ASSERT_EQUAL(1776411, hashed);
}

void test___onf_hash_rya_int_array___should_BigHashTest1(void)
{
  // I use signed ints, so I can handle kmers up to k length 15 without overflowing.
// ggggggggggggggg
  rya_int length = 15;
  rya_int ints[] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
  struct onf_rya_int_array* ary = onf_rya_int_array_new(15);
  ary->array = (rya_int * ) & ints;

  rya_int hashed = onf_hash_rya_int_array(ary);
  TEST_ASSERT_EQUAL(1073741823, hashed);
}

void test___onf_hash_rya_int_array___should_ReturnErrorOnOverflow(void)
{
  rya_int len    = 16; // should be big enough to overflow... 4 ** 16
  rya_int ints[] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};

  struct onf_rya_int_array* ary = onf_rya_int_array_new(len);
  ary->array = (rya_int * ) & ints;

  rya_int hashed = onf_hash_rya_int_array(ary);

  TEST_ASSERT_EQUAL(ONF_ERROR_INT32_T, hashed);
}


//////////////////////

void test___onf_kmer_count_array_new___should_ReturnNewKmerCountArray(void)
{
  size_t  size        = 2;
  size_t  output_size = 16;
  rya_int expected[]  = {
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
  };

  struct onf_rya_int_array* actual = onf_kmer_count_array_new(size);

  TEST_ASSERT_EQUAL(output_size, actual->length);
  TEST_ASSERT_EQUAL_INT_ARRAY(expected, actual->array, output_size);

  free(actual);
}

void test___onf_kmer_count_array_new___should_ReturnErrorIfSizeIsBad(void)
{
  size_t size = 0;

  struct onf_rya_int_array* actual = onf_kmer_count_array_new(size);

  TEST_ASSERT_EQUAL_PTR(ONF_ERROR_PTR, actual);

  free(actual);
}

//////////////////////

void test___onf_hash_lower_order_kmer___should_ReturnHashValOfLowerOrderKmer(void)
{
  rya_int how_much_lower = 1;

  struct onf_rya_int_array* kmer3 = onf_rya_int_array_new(3);
  kmer3->array[0] = 0;
  kmer3->array[1] = 1;
  kmer3->array[2] = 2;

  struct onf_rya_int_array* kmer2 = onf_rya_int_array_new(2);
  kmer2->array[0] = 0;
  kmer2->array[1] = 1;

  rya_int kmer3_hash = onf_hash_rya_int_array(kmer3);

  rya_int expected = onf_hash_rya_int_array(kmer2);

  rya_int actual = onf_hash_lower_order_kmer(kmer3_hash, how_much_lower);

  TEST_ASSERT_EQUAL_INT(expected, actual);
}

void test___onf_hash_lower_order_kmer___HandlesAnyLowerOrder(void)
{
  rya_int how_much_lower = 2;

  struct onf_rya_int_array* kmer3 = onf_rya_int_array_new(3);
  kmer3->array[0] = 0;
  kmer3->array[1] = 1;
  kmer3->array[2] = 2;

  struct onf_rya_int_array* kmer1 = onf_rya_int_array_new(1);
  kmer1->array[0] = 0;

  rya_int kmer3_hash = onf_hash_rya_int_array(kmer3);

  rya_int expected = onf_hash_rya_int_array(kmer1);

  rya_int actual = onf_hash_lower_order_kmer(kmer3_hash, how_much_lower);

  TEST_ASSERT_EQUAL_INT(expected, actual);
}

//////////////////////

void test___onf_count_kmers___should_ReturnKmerCountsForSequence(void)
{
  // kmers: ac, ct, tg, ga, ac
  // encoded kmers are {0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 1}
  char* seq = "actNgac"; // {0, 1, 2, 3, 0, 1 }

  size_t kmer_size = 2, seq_len = 7;

  rya_int counts[] = {
      0, 2, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1,
      1, 0, 0, 0,
  };

  struct onf_rya_int_array* actual = onf_count_kmers(seq, seq_len, kmer_size);

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

  rya_int num_kmers_in_seq     = seq_len - ksize + 1;
  rya_int total_possible_kmers = pow(4, ksize);

  struct onf_rya_int_array* actual_counts = onf_count_kmers(seq, seq_len, ksize);

  rya_int* expected_counts = calloc(total_possible_kmers, sizeof(int));
  // they are all the same kmer "a" * ksize except the last one
  expected_counts[0] = num_kmers_in_seq - 1;
  // this one is ("a" * ksize-1) + "c", which will hash to 1.
  expected_counts[1] = 1;

  TEST_ASSERT_EQUAL(total_possible_kmers, actual_counts->length);
  TEST_ASSERT_EQUAL_INT_ARRAY(expected_counts, actual_counts->array, total_possible_kmers);

  onf_rya_int_array_free(actual_counts);
  free(expected_counts);
}

void test___onf_count_kmers___should_ReturnErrorOnBadInput(void)
{
  struct onf_rya_int_array* actual = NULL;

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

//void print_rya_int_array(struct onf_rya_int_array* ary)
//{
//  printf("ary (%zu):", ary->length);
//  for (size_t i = 0; i < ary->length; ++i) {
////    printf(" %d", ary->array[i]);
//  }
//  putchar('\n');
//}
//
//void print_non_zero(struct onf_rya_int_array* ary, char* msg)
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

  struct onf_rya_int_array* counts6 = onf_count_kmers(seq, seq_len, 6);
  struct onf_rya_int_array* counts8 = onf_count_kmers(seq, seq_len, 8);
  struct onf_rya_int_array* counts9 = onf_count_kmers(seq, seq_len, 9);
  assert(counts6->length == pow(4, 6));
  assert(counts8->length == pow(4, 8));
  assert(counts9->length == pow(4, 9));

  onf_count_kmers2(seq, seq_len);

  // The actual function being tested.
  struct onf_rya_int_array** arrays = onf_count_kmers2(seq, seq_len);

  struct onf_rya_int_array* actual_counts9 = arrays[2];
  struct onf_rya_int_array* actual_counts8 = arrays[1];
  struct onf_rya_int_array* actual_counts6 = arrays[0];

  TEST_ASSERT_EQUAL(counts9->length, actual_counts9->length);
  TEST_ASSERT_EQUAL(counts8->length, actual_counts8->length);
  TEST_ASSERT_EQUAL(counts6->length, actual_counts6->length);

  TEST_ASSERT_EQUAL_INT_ARRAY(counts9->array, actual_counts9->array, counts9->length);
  TEST_ASSERT_EQUAL_INT_ARRAY(counts8->array, actual_counts8->array, counts8->length);
  TEST_ASSERT_EQUAL_INT_ARRAY(counts6->array, actual_counts6->array, counts6->length);

  onf_rya_int_array_free(counts6);
  onf_rya_int_array_free(counts8);
  onf_rya_int_array_free(counts9);

  for (rya_int i = 0; i < 2; ++i) {
    onf_rya_int_array_free(arrays[i]);
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