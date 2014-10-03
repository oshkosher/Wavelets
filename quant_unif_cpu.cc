#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dwt_cpu.h"

// Applies the threshold such that abs(values) <= threshold are 0
// Maps the remaining range of values to the +/- 0:(2^(bits-1))-1
// Overwrites data with the new values
void quant_unif_cpu(int len, float *data, int bits, float threshold, float maxVal)
{
	int count = len * len;
	// 2^(bits-1) - 1
	int base = (1 << (bits-1)) - 1;

	for (int idx = 0; idx < count; idx++ )
	{
		if (fabsf(data[idx]) <= threshold)
		{
			data[idx] = (float)0.0;
		}
		else
		{
			int sign=data[idx]/fabsf(data[idx]);
			float quantized = sign*ceil(base*((fabsf(data[idx])-threshold)/(maxVal-threshold)));
			// if (idx < 20) printf("data[%d] %f -> %f\n", idx, data[idx], quantized);
			data[idx] = quantized;
		}
	}
}
