#ifndef __QUANT_COUNT_H__
#define __QUANT_COUNT_H__

#include <cstdlib>
#include <cmath>
#include <cassert>
#include <vector>
#include <algorithm>

/**
   Generate a quantization based on bins containing equal numbers of
   input values.

   1. Sort the data by absolute value.
   2. Consider every datum with an absolute value less than equal to
      'thresh' as bin 0.
   3. Starting with the values greater than 'thresh', split them into
      (2^bits)-1 bins, each containing an equal number of values.
   4. The codebook value for the first bin (which contains 0) will
      be 0. The codebook value for each remaining bin will be the
      median value of all the datums in that bin.

   For example, for the given quantization:
     [0..1)  [1..5)  [5..10)  [10..inf)

   The bin boundaries would be [1, 5, 10].
   And a possible codebook would be [0, 2, 7, 20].

   If the number of elements below the threshold are at least the
   size of a bin, then the first bin threshold will be the given
   value 'thresh'. Otherwise, the first threshold will be larger than
   'thresh'.
*/

/**
   Given 'len' data points in data[] (not len*len), compute the bin
   boundaries and codebook to quantize the data into 2^bits bits.
*/
void quant_count_cpu(int len, float data[], int bits, float thresh,
		   std::vector<float> &boundaries,
		   std::vector<float> &codebook);

/** Like quant_pop_cpu, but the given data is already the sorted
    absolute values. */
void quant_count_cpu_sorted(int len, float sorted[], int bits, float thresh,
			  std::vector<float> &boundaries,
			  std::vector<float> &codebook);

#endif // __QUANT_COUNT_H__
