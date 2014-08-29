#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <math.h>
#include "dwt_cpu.h"

// Applies the threshold such that values <= threshold are 0
// Maps the remaining range of values to the values 0:(2^bits)-1
// Overwrites data with the new values
void dquant_log_cpu(int len, float *data, int bits, float threshold, float maxVal)
{
	int count = len * len;
	int base = pow(2,bits-1)-1;
     lmax = log2(maxVal/threshold);
	for (int idx = 0; idx < count; idx++ )
	{
		if (data[idx] <= threshold)
		{
			data[idx] = (float)0.0;
		}
		else
		{
			sign=data[idx]/abs(data[idx]);
			lnVal=(abs(data[idx])*lmax0/base;
			data[idx] = sign*threshold*(2^lnVal);
		}
	}
}
