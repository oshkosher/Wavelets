#ifndef __TEST_COMPRESS_CPU_H__
#define __TEST_COMPRESS_CPU_H__

#include "test_compress_common.h"
#include "wavelet.h"

// quantize the data
bool quantize(const CubeFloat &data, CubeInt &quantizedData,
              float maxAbsVal, WaveletCompressionParam &param,
              const float *nonzeroData, int nonzeroCount,
              float minValue, float maxValue,
              int *zeroBin = NULL);



#endif // __TEST_COMPRESS_CPU_H__
