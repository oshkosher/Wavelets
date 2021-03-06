#include <thrust/binary_search.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include "cucheck.h"
#include "optimize.h"
#include "test_compress_common.h"
#include "test_compress_gpu.h"
#include "cuda_timer.h"
#include "quant_gpu.h"
#include "nixtimer.h"
#include "histogram_gpu.h"
#include "huffman.h"


// return a value from 'sorted'
float OptimizationData::getSorted(int index) {
  if (index < 0 || index >= count()) {
    assert(index >= 0 && index < count());
    return 0;
  }

  float value;
  CUCHECK(cudaMemcpy(&value, sorted + index, sizeof(float),
                     cudaMemcpyDeviceToHost));

  return value;
}


// Look up the index of an item in 'sorted'
int OptimizationData::findSorted(float value) {
  thrust::device_ptr<const float> start(sorted), end(sorted+count());

  return thrust::lower_bound(start, end, value) - start;
}


__global__ void applyThresholdKernel(float *results, const float *input,
                                     int count, float thresholdValue) {
  for (int i = threadIdx.x + blockIdx.x*blockDim.x;
       i < count; i += blockDim.x*gridDim.x) {
    
    float f = input[i];
    if (fabsf(f) <= thresholdValue) f = 0;
    results[i] = f;
  }
}


/** Call this to try out compression parameters thresholdValue and binCount.
*/
bool testParameters(OptimizationData *o,
                    float thresholdValue, int binCount,
                    QuantizeAlgorithm quantAlg,
                    int *outputSizeBytes,
                    float *l1Error, float *l2Error, float *mse, float *pSNR,
                    float *relativeError) {

  WaveletCompressionParam param = o->transformedData->param;
  param.binCount = binCount;
  param.quantAlg = quantAlg;
  param.thresholdValue = thresholdValue;
  // param.originalSize = o->originalData->size;

  *outputSizeBytes = 0;

  int count = o->count();
  const float *inputData = o->transformedData->pointer(0,0,0);

  float *inverseWaveletInput_dev;
  CUCHECK(cudaMalloc((void**)&inverseWaveletInput_dev, count * sizeof(float)));

  // If binCount is nonpositive, don't quantize and dequantize.
  // Just apply the threshold and reverse the wavelet transform.
  if (binCount <= 0) {

    // apply threshold
    // copy the data to a new array, zeroing out values less than threshold

    CudaTimer timer("Apply threshold");
    timer.start();
    applyThresholdKernel<<<16,1024>>>
      (inverseWaveletInput_dev, (float*)o->transformedData->data_, count,
       thresholdValue);
    timer.end();
    if (!QUIET) {
      timer.sync();
      timer.print();
    }

  } else {  // binCount > 0

    // in one pass over the data:
    //   read one element
    //   quantize value
    //   increment appropriate histogram bin
    //   dequantize value
    //   write to new array
    // quantHistDequantKernel<<<16,1024>>>
    // (inverseWaveletInput_dev, (float*)transformedData.data_, count, 


    // but for comparison, first implement it with multiple passes

    int nonzeroIdx = o->findSorted(param.thresholdValue);
    const float *nonzeroData = o->sorted + nonzeroIdx;
    int nonzeroCount = count - nonzeroIdx;
    CudaTimer quantTimer("Quantize"), histTimer("Histogram"),
      dequantTimer("Dequantize"), histCopyTimer("Copy histogram data to CPU");
    int zeroBin;
    
    if (!quantizeGPU((int*)inverseWaveletInput_dev, inputData, count,
                     param, nonzeroData, nonzeroCount,
                     o->maxAbsVal, o->minVal, o->maxVal, quantTimer, &zeroBin))
      return false;

    /*
    printDeviceArray((int*)inverseWaveletInput_dev,
                     o->transformedData->size, "after quantize");
    */
    
    int *freqCounts;
    freqCounts = computeFrequenciesGPU
      (binCount, (int*)inverseWaveletInput_dev, count, zeroBin, false);

    // if the quantized data needs to be scanned for zero compression,
    // copy it to the CPU for processing
    int *quantizedData = NULL;
    if (o->doCompressZeros) {
      quantizedData = new int[count];
      CUCHECK(cudaMemcpy(quantizedData, inverseWaveletInput_dev,
                         sizeof(int) * count, cudaMemcpyDeviceToHost));
    } else {
      zeroBin = -1;
    }

    if (!dequantizeGPU(inverseWaveletInput_dev,
                       (const int*)inverseWaveletInput_dev, count, param))
      return false;

    // compute the huffman coding
    Huffman huff;
    initHuffman(huff, count, quantizedData, binCount, freqCounts, zeroBin);

    // just get the size of the encoded data; don't write it out
    *outputSizeBytes = huff.totalEncodedLengthBytes();

    if (!QUIET) {
      histCopyTimer.sync();
      quantTimer.print();
      // histTimer.print();
      // dequantTimer.print();
      // histCopyTimer.print();
      fflush(stdout);
    }
    
    delete[] freqCounts;
    if (quantizedData) delete[] quantizedData;
  }

  // compare to original
  ErrorAccumulator errAccum;
  errAccum.setMaxPossible(o->originalData->getMaxPossibleValue());

  float *tempData_dev;
  CUCHECK(cudaMalloc((void**)&tempData_dev, count * sizeof(float)));

  switch (o->originalData->datatype) {
  case WAVELET_DATA_UINT8:
    computeErrorRatesAfterDequantGPU
      (inverseWaveletInput_dev, o->transformedData->size, tempData_dev, param,
       (const unsigned char *)o->originalData->data_, o->originalData->datatype,
       errAccum);
       
    break;
  case WAVELET_DATA_INT32:
    computeErrorRatesAfterDequantGPU
      (inverseWaveletInput_dev, o->transformedData->size, tempData_dev, param,
       (const int *)o->originalData->data_, o->originalData->datatype,
       errAccum);
    break;
  case WAVELET_DATA_FLOAT32:
    computeErrorRatesAfterDequantGPU
      (inverseWaveletInput_dev, o->transformedData->size, tempData_dev, param,
       (const float *)o->originalData->data_, o->originalData->datatype,
       errAccum);
    break;
  default:
    break;
  }

  CUCHECK(cudaFree(tempData_dev));
  CUCHECK(cudaFree(inverseWaveletInput_dev));

  if (l1Error) *l1Error = errAccum.getL1Error();
  if (l1Error) *l2Error = errAccum.getL2Error();
  if (mse) *mse = errAccum.getMeanSquaredError();
  if (pSNR) *pSNR = errAccum.getPeakSignalToNoiseRatio();
  if (relativeError) *relativeError = errAccum.getRelativeError();

  return true;
}
