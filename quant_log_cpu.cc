#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include "quant_log_cpu.h"

float log2(float x)
{
return (log(x)/log(2.0));
}

// Applies the threshold such that values <= threshold are 0
// Maps the remaining range of values to the values 0:(2^bits)-1
// Overwrites data with the new values
float quant_log_cpu(int len, float *data, int bits, float threshold, float maxVal)
{
	int count = len * len;
	int base = (int)(pow(2.0,bits-1)-1);
    float lmax = (float)(log2((float)(maxVal/threshold)));
	for (int idx = 0; idx < count; idx++)
	{
		if (data[idx] <= threshold)
		{
			data[idx] = (float)0.0;
		}
		else
		{
			int sign=data[idx]/abs(data[idx]);
			float lnVal=log2((float(abs(data[idx])/threshold)));
			data[idx] = sign*ceil((base*lnVal)/lmax);
		}
	}
	
	return lmax;
}
