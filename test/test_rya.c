#include "unity.h"
#include "rya.h"

void setUp(void)
{
}

void tearDown(void)
{
}

void test___rya_file_exist___should_ReturnTrueIfFileExists(void)
{
  char* fname = __FILE__;
  TEST_ASSERT_EQUAL(rya_true, rya_file_exist(fname));
}

void test___rya_file_exist___should_ReturnFalseIfFileDoesntExist(void)
{
  char* fname = "arstoienarsotienarstoin.txt";
  TEST_ASSERT_EQUAL(rya_false, rya_file_exist(fname));
}

void test___rya_file_exist___should_ReturnErrorIfFnameIsBad(void)
{
  char* fname = NULL;
  TEST_ASSERT_EQUAL(rya_error, rya_file_exist(fname));
}
