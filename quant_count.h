#ifndef __QUANT_COUNT_H__
#define __QUANT_COUNT_H__

#include <cstdlib>
#include <cmath>
#include <cassert>
#include <vector>


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
