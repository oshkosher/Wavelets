#ifndef __DQUANT_LOG_CPU__
#define __DQUANT_LOG_CPU__

// Dequantizes the compressed data and outputs the dquantized data
// The dquantized data is allocated here and must later be deleted.
void dquant_log_cpu(int len, float *data, int bits, float threshold, float maxVal);

#endif // __DQUANT_LOG_CPU__
