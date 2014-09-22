#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <math.h>
#include "dwt_cpu.h"

// Applies the threshold such that values <= threshold are 0
// Maps the remaining range of values to the values 0:(2^bits)-1
// Overwrites data with the new values
void dquant_unif_cpu(int len, float *data, int bits, float threshold, float maxVal)
{
	int count = len * len;
	int base = pow(2.0,bits-1)-1;
	for (int idx = 0; idx < count; idx++ )
	{
		if (abs(data[idx]) <= threshold)
		{
			int sign=data[idx]/abs(data[idx]);
			data[idx] = sign*((abs(data[idx])*(maxVal/base)));			
		}
	}
}
