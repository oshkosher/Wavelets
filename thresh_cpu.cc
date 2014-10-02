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
   absolute value in the data.

   If saveSortedAbsData is non-NULL, then the caller wants a copy of the
   sorted array of absolute values. Set *saveSortedAbsData to point to
   the array, and do not deallocate it.
*/
float thresh_cpu(int len, float *data, float compRatio, float *maxVal, 
		 float **saveSortedAbsData)
{
	int count = len*len;	//Square matrix
	float *absData = new float[count];
	assert(absData);

	for (int idx = 0; idx < count; idx++ )
	{
		float afval = fabsf(data[idx]);
		absData[idx] = afval;
		// if (idx < 10) std::cout << "idx: " << idx << " data[idx]:" << data[idx] << " absData[idx]: " << afval << std::endl;
	}
	
	std::sort(absData, absData + count);

	int threshIdx = (int)(compRatio * count);
	float threshold = absData[threshIdx];
	if (threshold == 0) threshold = 1e-16;
	
	// std::cout << "Front:" << absData.front() << " Back:" << absData.back() << " threshold: " << threshold << std::endl;
	*maxVal = absData[count-1];   // get the largest value in the data
	
	// std::cout << "len:" << len << " maxVal:" << *maxVal << " threshIdx:" << threshIdx << " count:" << count << " threshold:" << threshold << std::endl;

	if (saveSortedAbsData) {
		*saveSortedAbsData = absData;
	} else {
		delete[] absData;
	}

	return threshold;
}
