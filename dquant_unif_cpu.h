#ifndef __QUANT_UNIF_CPU__
#define __QUANT_UNIF_CPU__

// Applies the threshold such that values <= threshold are 0
// Maps the remaining range of values to the values 0:(2^bits)-1
// Overwrites data with the new values
void quant_unif_cpu(int len, float *data, int bits, float threshold, float maxVal);

#endif // __QUANT_UNIF_CPU__
