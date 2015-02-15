#ifndef __HISTOGRAM_GPU_H__
#define __HISTOGRAM_GPU_H__

/**
   Compute frequency counts for values in an input array of integers.

   If onDevice is true, allocate the array on the device. The callers
   is responsible for freeing it with cudaFree().
   If onDevice is false, copy the data to the host. The callers
   is responsible for freeing it with delete[].

   binCount - number of bins

   data_dev - the input data, an array of integers
   
   count - length of data_dev array

   zeroBin - index of a bin that is likely to be very popular.
     The code optimizes for this bin by keeping its counter local to
     each thread, which greatly reduces contention.
*/
int *computeFrequenciesGPU(int binCount, const int *data_dev, int count,
                           int zeroBin, bool onDevice = true);


#endif // __HISTOGRAM_GPU_H__
