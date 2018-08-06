#include "unity.h"
#include "const.h"
#include "array.h"

void setUp(void)
{
}

void tearDown(void)
{
}

void test___onf_int_array_new___should_ReturnNewIntArray(void)
{
  size_t length = 3;
  struct onf_int_array* int_array = onf_int_array_new(length);

  int expected_ary[] = { 0, 0, 0 };

  TEST_ASSERT_EQUAL(length, int_array->length);
  TEST_ASSERT_EQUAL_INT_ARRAY(expected_ary, int_array->array, length);
}

void test___onf_int_array_new___should_ReturnErrorIfLenLessThanOne(void)
{
  size_t length = 0;
  struct onf_int_array* ary = onf_int_array_new(length);

  TEST_ASSERT_EQUAL_PTR(ONF_ERROR_PTR, ary);
}

////////////////////

void test___onf_int_array_bad___should_ReturnTrueIfArrayIsBad(void)
{
  size_t length = 10;
  struct onf_int_array* ary = onf_int_array_new(length);

  TEST_ASSERT_FALSE(onf_int_array_bad(ary));

  ary->array = NULL;
  TEST_ASSERT_TRUE(onf_int_array_bad(ary));

  ary = NULL;
  TEST_ASSERT_TRUE(onf_int_array_bad(ary));
}