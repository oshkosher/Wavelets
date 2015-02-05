#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <cfloat>
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

  // compute the absolute values of every element, and at the same time
  // find the largest (most positive) and smallest (most negative) values
  for (int idx = 0; idx < count; idx++) {
    // printf("[%d] = %.8g %g %g\n", idx, data[idx], *minVal, *maxVal);
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
  if (threshold == 0) threshold = MIN_THRESHOLD_VALUE;

  *nonzeroCount = count - threshIdx - 1;

  if (saveSortedAbsData) {
    *saveSortedAbsData = absData;
  } else {
    delete[] absData;
  }

  return threshold;
}

class Traverse {
public:
  float min, max, *absData;
  int count;

  Traverse(float *a) : absData(a) {
    count = 0;
    min = FLT_MAX;
    max = FLT_MIN;
  }

  void visit(float x) {
    // printf("[%d] = %.8g %f %f\n", count, x, min, max);
    absData[count++] = fabsf(x);
    if (x < min) min = x;
    if (x > max) max = x;
  }
};


float thresh_cpu(CubeFloat *data, float compRatio,
		 int *nonzeroCount, float *maxVal, float *minVal,
		 float **saveSortedAbsData) {
		 
  float *absData = new float[data->count()];
  assert(absData);

  // compute the absolute values of every element, and at the same time
  // find the largest (most positive) and smallest (most negative) values
  Traverse traverser(absData);
  data->visit<Traverse>(traverser);
  *minVal = traverser.min;
  *maxVal = traverser.max;
	
  std::sort(absData, absData + data->count());

  int threshIdx = (int)(compRatio * data->count());
  float threshold = absData[threshIdx];
  if (threshold == 0) threshold = 1e-16;

  *nonzeroCount = data->count() - threshIdx - 1;

  if (saveSortedAbsData) {
    *saveSortedAbsData = absData;
  } else {
    delete[] absData;
  }

  return threshold;
}
