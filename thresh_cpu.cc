#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <math.h>
#include <vector>
#include <algorithm>    // std::sort
#include "thresh_cpu.h"

// Allocates absData on the heap and loads it with abs(data)
// Sorts absData.  Finds the index corresponding to the compRatio, and 
// reads that element to get the threshold. Sets min and max values that 
// remain above threshold  
float thresh_cpu(int len, float *data, float compRatio, float *maxVal)
{
	int count = len*len;	//Square matrix
	float *absData = new float[count];
	for (int idx = 0; idx < count; idx++ )
	{
		absData[idx] = abs(data[idx]);
	}
	
	int threshIdx = floor(compRatio * count );
	std::vector<float> dvec(absData, absData + count);
	std::sort(dvec.begin(),dvec.end());
	*maxVal = *(dvec.rbegin());  // get the largest value in the data
	float threshold = dvec[threshIdx];

	delete[] absData;
	return threshold;
}
