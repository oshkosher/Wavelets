#ifndef __TEST_COMPRESS_GPU_H__
#define __TEST_COMPRESS_GPU_H__

#include "wavelet.h"
#include "test_compress_common.h"

// return the number of mismatches
int compareArrays(const float *a, const float *b, int count);
int compareArrays(const int *a, const int *b, int count);

// debug output
void printArray(const float *array, int width, int height, int depth,
                const char *name = NULL);
void printArray(const int *array, int width, int height, int depth,
                const char *name = NULL);
void printDeviceArray(const float *array_dev, int width, int height, int depth,
                      const char *name = NULL);
void printDeviceArray(const float *array_dev, scu_wavelet::int3 size,
                      const char *name = NULL);
void printDeviceArray(const int *array_dev, int width, int height, int depth,
                      const char *name = NULL);
void printDeviceArray(const int *array_dev, scu_wavelet::int3 size,
                      const char *name = NULL);

void computeErrorRatesAfterDequantGPU
(float *data_dev, scu_wavelet::int3 size, float *tempData_dev,
 const WaveletCompressionParam &param,
 const unsigned char *inputData_dev, ErrorAccumulator &errAccum);


#endif // __TEST_COMPRESS_GPU_H__
