#ifndef __QUANT_LOG_CPU__
#define __QUANT_LOG_CPU__

// Define log2 on systems that don't have it
float quant_log2(float x);

// Applies the threshold such that values <= threshold are 0
// Maps the remaining range of values to the values 0:(2^bits)-1
// Overwrites data with the new values
float  quant_log_cpu(int len, float *data, int bits, float threshold, float maxVal);


#endif // __QUANT_LOG_CPU__
