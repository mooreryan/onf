#include <assert.h>
#include <math.h>
#include <rya_file.h>

#include "onf.h"
#include "onf_math.h"
#include "array.h"
#include "tommyarray.h"

double score(struct onf_rya_int_array** host_counts,
             struct onf_rya_int_array** virus_counts,
             double pseudo_count)
{
  assert(pseudo_count > 0);
  // TODO assert all the rest of the stuff too.

  // We need the 8mer counts and the 9mer counts.
  //  int sizes[]   = {6, 8, 9};

  struct onf_rya_int_array *host_counts8, *host_counts9, *virus_counts8, *virus_counts9;

  host_counts8 = host_counts[1];
  host_counts9 = host_counts[2];

  virus_counts8 = virus_counts[1];
  virus_counts9 = virus_counts[2];

  // v9mer: read viral 9mer hash value
  // h9mer: read host 9mer hash value
  size_t v9mer = 0, v8mer = 0, h9mer = 0, h8mer = 0;

  rya_int h9mer_count = 0, h8mer_count = 0, v8mer_count = 0;

  double sum = 0.0, val = 0.0;

  int non_zero_v8mers = 0;

  for (size_t v9mer = 0; v9mer < virus_counts9->length; ++v9mer) {
    v8mer = onf_hash_lower_order_kmer(v9mer, 1);

    assert(v9mer < host_counts9->length);
    assert(v8mer < host_counts8->length);

    // We only want to calculate score for kmers actually existing in the virus.
    v8mer_count = virus_counts8->array[v8mer];
    if (v8mer_count > 0) {
      ++non_zero_v8mers;
      h9mer_count = host_counts9->array[v9mer];
      h8mer_count = host_counts8->array[v8mer];

      val = (h9mer_count + 0.25 * pseudo_count) / (h8mer_count + pseudo_count);

      sum += log(val);
    }
  }

  // TODO should this be 8 or 9?

  // Return the mean
  return sum / non_zero_v8mers;
}

