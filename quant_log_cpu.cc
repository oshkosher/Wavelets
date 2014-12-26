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

  int count = len * len;
  int base = (1 << (bits-1)) - 1;
  float lmax = logf(maxVal/threshold);
  for (int idx = 0; idx < count; idx++) {
    float before = data[idx];
    float absdata = fabsf(data[idx]);

    if (absdata <= threshold) {
      data[idx] = 0.0f;
    } else {
      int sign=data[idx]/absdata;
	
      float lnVal=logf(absdata/threshold);
      data[idx] = sign*ceil((base*lnVal)/lmax);

    }
    printf("%f\t%.0f\n", before, data[idx]);
  }

  return lmax;
}
