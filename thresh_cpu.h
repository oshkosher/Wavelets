#ifndef __THRESH_CPU__
#define __THRESH_CPU__

// Allocates absData on the heap and loads it with abs(data)
// Sorts absData.  Finds the index corresponding to the compRatio, and 
// reads that element to get the threshold. Sets min and max values that 
// remain above threshold
float thresh_cpu(int len, float *data, float compRatio, float *maxVal, float *minVal) ;


#endif // __THRESH_CPU__
