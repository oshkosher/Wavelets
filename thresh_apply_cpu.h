#ifndef __THRESH_APPLY_CPU__
#define __THRESH_APPLY_CPU__

// Applies the threshold such that values <= threshold are 0
// Maps the remaining range of values to the values 0:(2^bits)-1
// Overwrites data with the new values
void thresh_apply_cpu(int len, float *data, int bits, float threshold, float maxVal);

#endif // __THRESH_APPLY_CPU__
