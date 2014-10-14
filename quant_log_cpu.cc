#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include "quant.h"
#include "quant_log_cpu.h"

// Applies the threshold such that values <= threshold are 0
// Maps the remaining range of values to the values 0:(2^bits)-1
// Overwrites data with the new values
float quant_log_cpu(int len, float *data, int bits, float threshold, float maxVal)
{
	int displayCount = 0;
	int count = len * len;
	int base = (int)(pow(2.0,bits-1)-1);
    float lmax = (float)(quant_log2((float)(maxVal/threshold)));
	for (int idx = 0; idx < count; idx++)
	{
		if (fabsf(data[idx]) <= threshold)
		{
			data[idx] = (float)0.0;
		}
		else
		{
			int sign=data[idx]/fabsf(data[idx]);
			
			float lnVal=quant_log2((float)(fabsf(data[idx]/threshold)));
			data[idx] = sign*ceil((base*lnVal)/lmax);
			if ( displayCount > 0 )
			{
				displayCount--;
				std::cout << "Idx: " << idx << " Value: " << data[idx] << std::endl;
				std::cout << "    sign: " << sign << std::endl;
				std::cout << "    lnVal: " << lnVal << std::endl;
				std::cout << "    threshold: " << threshold << " maxVal: " << maxVal << std::endl;
			}
		}
	}
	
	return lmax;
}
