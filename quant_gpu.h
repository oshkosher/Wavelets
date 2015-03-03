#ifndef __QUANT_GPU_H__
#define __QUANT_GPU_H__

#include "cuda_timer.h"
#include "wavelet.h"

bool quantizeGPU(int *outputData, const float *inputData, int count,
                 WaveletCompressionParam &param,
                 const float *nonzeroData_dev, int nonzeroCount,
                 float maxAbsVal, float minValue, float maxValue,
                 CudaTimer &quantizeTimer, int *zeroBin = NULL);

bool dequantizeGPU(float *result_dev, const int *input_dev,
                   int count, const WaveletCompressionParam &param);

void computeLloydQuantizationGPU(const float *inputData, int count,
                                 int binCount, float minVal, float maxVal,
                                 float thresholdValue,
                                 std::vector<float> &quantBinBoundaries,
                                 std::vector<float> &quantBinValues);

void __global__ quantUniformKernel
(int *output, const float *input, int count,
 int binCount, float threshold, float max);

void __global__ quantLogKernel
(int *output, const float *input, int count,
 int binCount, float threshold, float max);

void __global__ quantCodebookKernel
(int *output, const float *input, int count, 
 const float * __restrict boundaries, int boundaryCount);

void quantCodebookGPU
(int *output, const float *input, int dataCount,
 const vector<float> &boundaries, int binCount);

void __global__ dequantUniformKernel
(float *output, const int *input, int count,
 int binCount, float threshold, float max);

void __global__ dequantLogKernel
(float *data, const int *input, int count,
 int binCount, float threshold, float max);

void __global__ dequantCodebookKernel
(float *output, const int *input, int count,
 const float *binValues, int binCount);

void dequantCodebookGPU
(float *output, const int *input, int count,
 const vector<float> &binValues);


#endif // __QUANT_GPU_H__
