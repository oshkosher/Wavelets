#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <math.h>
#include "dwt_cpu.h"

// Allocates absData on the heap and loads it with abs(data)
// Sorts absData.  Finds the index corresponding to the compRatio, and 
// reads that element to get the threshold. Sets min and max values that 
// remain above threshold
float threshold thresh_cpu(int len, float *data, float compRatio, float *maxVal, float *minVal)
{
	count = len*len;	//Square matrix
	absData = new float[count];
	for (int idx = 0; idx < count; idx++ )
	{
		absData[idx] = abs(data[idx]);
	}
	
	int threshIdx = floor(compRatio * count );
	std::vector<float> dvec(absData, absData + count);
	std::sort(dvec.begin(),dvec.end());
	*maxVal = dvec.rbegin();  // get the largest value in the absData
	threshold = dvec[threshIdx];
	minVal = threshold + 1e-6;  // Like an eps, incase it falls through.
	for (int cidx = threshIdx+1; cidx < count; cidx++ )
	{
		if (dvec[cidx] == threshold ) continue;   // Skip past any entries equal to threshold
		*minVal = dvec[cidx];
		break;
	}

	delete[] absData;
	return threshold;
}
