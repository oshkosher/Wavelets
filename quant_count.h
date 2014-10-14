#ifndef __QUANT_COUNT_H__
#define __QUANT_COUNT_H__

#include <cstdlib>
#include <cmath>
#include <cassert>
#include <vector>

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
void quant_count_init_cpu
(int len, const float data[], int bits, float thresh,
 std::vector<float> &boundaries,
 std::vector<float> &codebook);

/** Like quant_count_init, but the given data is already the sorted
    absolute values. */
void quant_count_init_sorted_cpu
(int len, const float sorted[], int bits, float thresh,
 std::vector<float> &boundaries,
 std::vector<float> &codebook);

/** Apply the bin boundaries to a value; return the quantized value. */
int quant_boundaries(const std::vector<float> &boundaries, float value);

/** Apply the bin boundaries to an array of values */
void quant_boundaries_array(const std::vector<float> &codebook,
			    int len, float *data);


/** Given a value in the range [-(codebook.size()-1) .. (codebook.size()-1)],
    return the value in the codebook and return it. Note that if the input
    value is negative, the codebook value needs to be inverted. */
float dequant_codebook(const std::vector<float> &codebook, int value);

void dequant_codebook_array(const std::vector<float> &codebook,
			    int len, float *data);

#endif // __QUANT_COUNT_H__
