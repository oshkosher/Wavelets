#include "quant.h"


/**
   Generate a quantization based on bins containing equal numbers of
   input values.

   1. Sort the data
   2. Consider every datum with an absolute value less than equal to
      'thresh' as a bin with a codebook value of 0.
   3. Allocate the remaining bins proportionally between the values
      left of 0 and right of 0.
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
void QuantCodebook::initCountBins(int count, float *data,
				  int bits, float thresh) {

  boundaries.clear();
  codebook.clear();
  
  int binCount = 1 << bits;
  float *sorted = new float[count];
  memcpy(sorted, data, count * sizeof(float));
  std::sort(sorted, sorted + count);

  // index of the first entry that is <= the threshold
  int negThreshIdx = std::lower_bound(sorted, sorted+count, -thresh)
    - sorted;

  // index of the first entry that is > the threshold
  int posThreshIdx = std::upper_bound(sorted, sorted+count, thresh)
    - sorted;

  // number of values that are negative or positive
  int negCount = negThreshIdx;
  int posCount = count - posThreshIdx;

  // number of bins representing negative and positive numbers
  int negBinCount = (float)negCount / (negCount+posCount) * binCount;
  int posBinCount = binCount - negBinCount - 1;

  // number of values in each bin
  float binPop = (float) negCount / negBinCount;

  // assign the negative bins
  float binLeftIdx, binRightIdxFrac;
  int binMidIdx, binRightIdx;

  // use the negative extrema as the first codebook entry
  boundaries.push_back(sorted[(int)binPop]);
  codebook.push_back(sorted[0]);
  binLeftIdx = binPop;

  for (int i=1; i < negBinCount-1; i++) {
    binRightIdxFrac = binLeftIdx + binPop;
    binRightIdx = (int) binRightIdxFrac;
    binMidIdx = (int) (binLeftIdx + .5f*binPop);

    boundaries.push_back(sorted[binRightIdx]);

    // use the median of each bin as the codebook entry
    codebook.push_back(sorted[binMidIdx]);

    binLeftIdx = binRightIdxFrac;
  }

  // last negative bin
  boundaries.push_back(-thresh);
  binMidIdx = (int) ((negThreshIdx + binLeftIdx) / 2);
  codebook.push_back(sorted[binMidIdx]);

  // zero bin
  boundaries.push_back(thresh);
  codebook.push_back(0);

  // positive bins
  binPop = (float) posCount / posBinCount;
  binLeftIdx = posThreshIdx;
  for (int i=0; i < posBinCount-1; i++) {
    binRightIdxFrac = binLeftIdx + binPop;
    binRightIdx = (int) binRightIdxFrac;
    binMidIdx = (int) (binLeftIdx + .5f*binPop);

    boundaries.push_back(sorted[binRightIdx]);

    // use the median of each bin as the codebook entry
    codebook.push_back(sorted[binMidIdx]);

    binLeftIdx = binRightIdxFrac;
  }
    
  // use the positive extrema as the last codebook entry
  codebook.push_back(sorted[count-1]);

  // and give it its own bin
  lastBoundary = .01f * sorted[count-2] + .99f * sorted[count-1];
  boundaries[boundaries.size()-1] = lastBoundary;

  delete[] sorted;
}


void QuantCodebook::printCodebook() {
  printf("Codebook:\n");

  for (size_t i = 0; ; i++) {
    printf("Bin %3d. %f\n", (int)i, codebook[i]);
    if (i == boundaries.size()) break;
    printf("  %f\n", boundaries[i]);
  }
}
