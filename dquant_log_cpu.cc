#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include "quant.h"
#include "dquant_log_cpu.h"


// Applies the threshold such that values <= threshold are 0
// Maps the remaining range of values to the values 0:(2^bits)-1
// Overwrites data with the new values
void dquant_log_cpu(int len, float *data, int bits, float threshold, float maxVal)
{
	int count = len * len;
	int base = (int)(pow(2.0,bits-1)-1);
        float lmax = (float)(QuantFunction::quant_log2((float)(maxVal/threshold)));
	
	for (int idx = 0; idx < count; idx++ )
	{
		if (data[idx] != 0)
		{
			int sign=data[idx]/fabsf(data[idx]);
			float lnVal=fabsf(data[idx]*(lmax/base));
			data[idx] = sign*threshold*(float)(pow((float)(2.0),lnVal));
		}
	}
}
