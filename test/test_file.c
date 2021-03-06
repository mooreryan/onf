#include "unity.h"
#include "file.h"
#include "tommyarray.h"
#include "const.h"
#include "rya.h"
#include "rlib.h"

#include <assert.h>

void setUp(void)
{
}

void tearDown(void)
{
}

void test___onf_file_files_in_dir___should_GiveNamesOfTheFiles(void)
{

  char* this_fname = __FILE__;

  rstring* this_file = rstring_new(this_fname);
  assert(this_file);

  rstring* dirname = rfile_dirname(this_file);
  assert(dirname);

  rstring* path = rstring_format("%s/test_files/apple", rstring_data(dirname));
  assert(path);

  tommy_array* expected_ary = malloc(sizeof(tommy_array));
  tommy_array_init(expected_ary);

  rstring* fname = NULL;

  fname = rstring_format("%s/good.fa", rstring_data(path));
  assert(fname);
  tommy_array_insert(expected_ary, fname);

  fname = rstring_format("%s/is.txt", rstring_data(path));
  assert(fname);
  tommy_array_insert(expected_ary, fname);

  fname = rstring_format("%s/pie/amelia.fa", rstring_data(path));
  assert(fname);
  tommy_array_insert(expected_ary, fname);

  fname = rstring_format("%s/pie/ryan.txt", rstring_data(path));
  assert(fname);
  tommy_array_insert(expected_ary, fname);



  tommy_array* actual = onf_file_files_in_dir(rstring_data(path));
  assert(actual);

  TEST_ASSERT_EQUAL(4, tommy_array_size(actual));

  for (int i = 0; i < 4; ++i) {
    TEST_ASSERT_EQUAL_STRING(rstring_data(tommy_array_get(expected_ary, i)),
                             tommy_array_get(actual, i));
  }

  for (int i = 0; i < 4; ++i) {
    free(tommy_array_get(expected_ary, i));
  }

  rstring_free(this_file);
  rstring_free(dirname);
  rstring_free(path);
  tommy_array_done(expected_ary);
  free(expected_ary);
}

void test___onf_file_files_in_dir___should_ReturnErrorOnErrors(void)
{
  char* path = NULL;

  TEST_ASSERT_EQUAL_PTR(ONF_ERROR_PTR, onf_file_files_in_dir(path));

  path = "arstoienarstoienarsotien"; // this doesn't exist.

  TEST_ASSERT_EQUAL_PTR(ONF_ERROR_PTR, onf_file_files_in_dir(path));
}