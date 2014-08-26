#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <math.h>
#include "dwt_cpu.h"

// Applies the threshold such that values <= threshold are 0
// Maps the remaining range of values to the values 0:(2^bits)-1
// Overwrites data with the new values
void quant_unif_cpu(int len, float *data, int bits, float threshold, float maxVal, float minVal)
{
	int count = len * len;
	int base = pow(2,bits)-1;
	float factor = base/(maxVal-minVal);
	
	std::vector<float> dvec(data, data + count);
	std::sort(dvec.begin(),dvec.end());
	threshold = dvec[threshIdx];
	
	for (int idx = 0; idx < count; idx++ )
	{
		if (data[idx] <= threshold)
		{
			data[idx] = (float)0.0;
		}
		else
		{
			data[idx] = ((data[idx]-minVal) * factor);
		}
	}
}
