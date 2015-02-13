#ifndef __HISTOGRAM_GPU_H__
#define __HISTOGRAM_GPU_H__

/**
   Compute frequency counts for values in an input array of integers.

   freqCounts_dev - an array of length 'binCount' will be allocated on
     the GPU, and the address will be assigned to this variable.
     The results will be in this array.

   binCount - number of bins

   data_dev - the input data, an array of integers
   
   count - length of data_dev array

   zeroBin - index of a bin that is likely to be very popular.
     The code optimizes for this bin by keeping its counter local to
     each thread, which greatly reduces contention.
*/
void computeFrequenciesGPU(int *&freqCounts_dev, int binCount,
                           const int *data_dev, int count, int zeroBin);


#endif // __HISTOGRAM_GPU_H__
