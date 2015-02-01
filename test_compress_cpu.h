#ifndef __TEST_COMPRESS_CPU_H__
#define __TEST_COMPRESS_CPU_H__

#include "test_compress_common.h"
#include "wavelet.h"

// quantize the data
bool quantize(const CubeFloat &data, CubeInt &quantizedData,
              float maxAbsVal, WaveletCompressionParam &param,
              const float *nonzeroData, int nonzeroCount,
              float minValue, float maxValue,
              float *quantErrorOut = NULL);

bool dequantize(const CubeInt &quantizedData, CubeFloat &data,
                const WaveletCompressionParam &param);

// this will modify restoredData in place
void computeErrorRatesAfterDequant
(CubeFloat *restoredData,
 const WaveletCompressionParam &param,
 const Cube *inputData,
 ErrorAccumulator *errAccum);



#endif // __TEST_COMPRESS_CPU_H__
