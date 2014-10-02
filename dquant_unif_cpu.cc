#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <math.h>
#include "dwt_cpu.h"

/**
  Maps the values 0:(2^bits)-1 to 0:maxVal.
  If the value is negative, the result will be negative.
  Overwrites data with the new values.
*/
void dquant_unif_cpu(int len, float *data, int bits, float threshold, float maxVal)
{
	int count = len * len;
	// 2^(bits-1) - 1
	int base = (1 << (bits-1)) - 1;
	for (int idx = 0; idx < count; idx++ )
	{
		float dequant = data[idx] * maxVal / base;
		// if (idx < 20) printf("data[%d] %f -> %f\n", idx, data[idx], dequant);
		data[idx] = dequant;
	}
}
