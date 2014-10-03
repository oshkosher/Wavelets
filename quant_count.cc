#include "quant_count.h"
#include <algorithm>


/**
   Given 'len' data points in data[] (not len*len), compute the bin
   boundaries and codebook to quantize the data into 2^bits bits.
*/
void quant_count_init_cpu
(int len, float data[], int bits, float thresh,
 std::vector<float> &boundaries,
 std::vector<float> &codebook) {
  
  // sort the data by absolute value
  float *sorted = new float[len];
  assert(sorted);
  for (int i=0; i < len; i++)
    sorted[i] = fabsf(data[i]);

  std::sort(sorted, sorted+len);

  quant_count_init_sorted_cpu(len, sorted, bits, thresh, boundaries, codebook);

  delete[] sorted;  
}

/** Like quant_count_cpu, but the given data is already the sorted
    absolute values. */
void quant_count_init_sorted_cpu
(int len, float sorted[], int bits, float thresh,
 std::vector<float> &boundaries,
 std::vector<float> &codebook) {

  assert(bits >= 1 && bits <= 32);

  int binCount = 1 << bits;
  // round down the bin size, so the last bin will be the largest.
  // for example, with 5 elements and 4 bins, the first three will be
  // of size 1, the last will contain 2.
  int binSize = len / binCount;

  boundaries.resize(binCount-1);
  codebook.resize(binCount);
  codebook[0] = 0;

  // is the flush-to-zero threshold smaller than the first bin?
  float *threshPtr = std::upper_bound(sorted, sorted+len, thresh);
  int threshPos = threshPtr - sorted;

  // index in sorted[] of the last value in the previous bin
  int binEnd;

  if (threshPos < binSize) {
    // the regular size of the first bin is enough to include everything
    // below the threshold value
    binEnd = binSize;
    boundaries[0] = sorted[binEnd];
  } else {
    // There are more values below the threshold than a normal-sized bin.
    // Make them their own bin, and resize the remainder accordingly.
    boundaries[0] = sorted[threshPos];
    binEnd = threshPos;
    binSize = (len - threshPos) / (binCount-1);
  }

  // There will always be 2^bits bins, and (2^bits)-1 boundaries.
  // The first bin is assumed to start at 0, so its lower bound is omitted.
  // Each boundary defines the lower bound of the bin.
  // Example: 4 bins the these boundaries:
  //   0......3.......5......7......
  // 2.99 will go in the first bin, but 3.00 will go in the second.

  for (int i=1; i < binCount-1; i++) {
    int binStart = binEnd;
    binEnd += binSize;

    // If there's a value that's so common it spans multiple bins, 
    // find the last copy of that value and use the value right after
    // it as the next boundary.
    if (binEnd < len && sorted[binEnd] == boundaries[i-1]) {
      float *next = std::upper_bound(sorted+binEnd, sorted+len, boundaries[i-1]);
      binEnd = next - sorted;
      binSize = (len - binEnd) / (binCount - i);
    }
    
    if (binEnd >= len) {
      boundaries[i] = codebook[i] = sorted[len-1];
      binSize = 0;
    } else {
      boundaries[i] = sorted[binEnd];
      codebook[i] = sorted[(binStart+binEnd)/2];
    }
  }
  codebook[binCount-1] = sorted[(binEnd + len)/2];

}

/** Apply the bin boundaries to a value; return the quantized value. */
int quant_boundaries(const std::vector<float> &boundaries, float value) {

  int sign = 1;
  if (value < 0) {
    sign = -1;
    value = -value;
  }

  // if the value is >= the end of the last bin, return the last bin
  int lastPos = (int)boundaries.size()-1;
  float lastBound = boundaries[lastPos];
  if (value >= lastBound)
    return sign * (lastPos+1);

  // Find the first boundary that is greater than value.
  // For example, if the boundaries are:
  // boundaries[0] = 3
  // boundaries[1] = 5
  // boundaries[2] = 10
  // Given .5, it returns 0 because 3 > .5
  // Given 5, it returns 2, because 10 > 5
  // Given 100, it return 3, because no entry is > 100
  
  std::vector<float>::const_iterator pos = 
    std::upper_bound(boundaries.begin(), boundaries.end(), value);
  int bin = (int) (pos - boundaries.begin());

  return sign * bin;
}
  

  

/** Apply the codebook to an array of values */
void quant_boundaries_array(const std::vector<float> &boundaries,
			    int len, float *data) {
  for (int i=0; i < len; i++) {
    data[i] = (float) quant_boundaries(boundaries, data[i]);
  }
}



/** Given a value in the range [-(codebook.size()-1) .. (codebook.size()-1)],
    return the value in the codebook and return it. Note that if the input
    value is negative, the codebook value needs to be inverted. */
float dequant_codebook(const std::vector<float> &codebook, int value) {
  assert(abs(value) < (int)codebook.size());

  if (value < 0) {
    return -codebook[-value];
  } else {
    return codebook[value];
  }
}



void dequant_codebook_array(const std::vector<float> &codebook,
			    int len, float *data) {
  for (int i=0; i < len; i++) {
    data[i] = dequant_codebook(codebook, data[i]);
  }
}


