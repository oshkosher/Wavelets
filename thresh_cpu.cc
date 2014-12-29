#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>    // std::sort
#include "thresh_cpu.h"

/**
   Allocates absData on the heap and loads it with abs(data)
   Sorts absData.  Finds the index corresponding to the compRatio, and 
   reads that element to get the threshold. Sets *maxVal to the largest
   value in the data and *minVal to the smallest (usually negative).

   If saveSortedAbsData is non-NULL, then the caller wants a copy of the
   sorted array of absolute values. Set *saveSortedAbsData to point to
   the array, and do not deallocate it.
*/
float thresh_cpu(int count, float *data, float compRatio,
		 int *nonzeroCount, float *maxVal, float *minVal,
		 float **saveSortedAbsData) {
		 
  float *absData = new float[count];
  assert(absData);

  *maxVal = *minVal = data[0];

  for (int idx = 0; idx < count; idx++) {
    absData[idx] = fabsf(data[idx]);
    if (data[idx] < *minVal) {
      *minVal = data[idx];
    } else if (data[idx] > *maxVal) {
      *maxVal = data[idx];
    }
  }
	
  std::sort(absData, absData + count);

  int threshIdx = (int)(compRatio * count);
  float threshold = absData[threshIdx];
  if (threshold == 0) threshold = 1e-16;

  *nonzeroCount = count - threshIdx - 1;

  if (saveSortedAbsData) {
    *saveSortedAbsData = absData;
  } else {
    delete[] absData;
  }

  return threshold;
}
