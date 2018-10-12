#ifndef __TEST_COMPRESS_GPU_H__
#define __TEST_COMPRESS_GPU_H__

#include "wavelet.h"
#include "test_compress_common.h"

// return the number of mismatches
int compareArrays(const float *a, const float *b, int count);
int compareArrays(const int *a, const int *b, int count);

template <class T>
void computeErrorRatesAfterDequantGPU
(float *data_dev, scu_wavelet::int3 size, float *tempData_dev,
 const WaveletCompressionParam &param,
 const T *inputData_dev, WaveletDataType inputType, ErrorAccumulator &errAccum);

#endif // __TEST_COMPRESS_GPU_H__
