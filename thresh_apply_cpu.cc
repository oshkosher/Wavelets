#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dwt_cpu.h"

// Applies the threshold such that abs(values) <= threshold are 0
// Maps the remaining range of values to the +/- 0:(2^(bits-1))-1
// Overwrites data with the new values
void thresh_apply_cpu(int len, float *data, int bits, float threshold, float maxVal)
{
	int count = len * len;
	int base = (int)(pow(2.0,bits-1)-1);

	for (int idx = 0; idx < count; idx++ )
	{
		if (abs(data[idx]) <= threshold)
		{
			data[idx] = (float)0.0;
		}
	}
}
