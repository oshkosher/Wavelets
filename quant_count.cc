#include "quant_count.h"


/**
   Given 'len' data points in data[] (not len*len), compute the bin
   boundaries and codebook to quantize the data into 2^bits bits.
*/
void quant_count_cpu(int len, float data[], int bits, float thresh,
		     std::vector<float> &boundaries,
		     std::vector<float> &codebook) {
  
  // sort the data by absolute value
  float *sorted = new float[len];
  assert(sorted);
  for (int i=0; i < len; i++)
    sorted[i] = fabsf(data[i]);

  std::sort(sorted, sorted+len);

  quant_count_cpu_sorted(len, sorted, bits, thresh, boundaries, codebook);

  delete[] sorted;  
}

/** Like quant_count_cpu, but the given data is already the sorted
    absolute values. */
void quant_count_cpu_sorted(int len, float sorted[], int bits, float thresh,
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
  int binEnd;
  if (threshPos < binSize) {
    // the regular size of the first bin is enough to include everything
    // below the threshold value
    binEnd = binSize;
    boundaries[0] = sorted[binEnd];
  } else {
    // there are more values below the threshold than a normal-sized bin
    // make them their own bin, and resize the remainder accordingly
    boundaries[0] = sorted[threshPos];
    binEnd = threshPos;
    binSize = (len - threshPos) / (binCount-1);
  }

  for (int i=1; i < binCount-1; i++) {
    int binStart = binEnd;
    binEnd += binSize;
    // if it's a really common value, find the next entry that's larger
    if (sorted[binEnd] == boundaries[i-1]) {
      float *next = std::upper_bound(sorted+binEnd, sorted+len, boundaries[i-1]);
      binEnd = next - sorted;
      binSize = (len - binEnd) / (binCount - i);
    }

    boundaries[i] = sorted[binEnd];
    codebook[i] = sorted[(binStart+binEnd)/2];
  }
  codebook[binCount-1] = sorted[(binEnd + len)/2];

}
