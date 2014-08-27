#ifndef __DQUANT_UNIF_CPU__
#define __DQUANT_UNIF_CPU__

// Dequantizes the compressed data and returns the new array.
// The output array is allocated here and must later be deleted.
float * quant_unif_cpu(int len, float *data, int bits, float threshold, float maxVal);

#endif // __DQUANT_UNIF_CPU__
