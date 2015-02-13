#include <thrust/binary_search.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include "cucheck.h"
#include "optimize.h"
#include "test_compress_common.h"

// from test_compress_gpu.cu
void computeErrorRatesAfterDequantGPU
(float *data_dev, scu_wavelet::int3 size, float *tempData_dev,
 const WaveletCompressionParam &param,
 const unsigned char *inputData_dev, ErrorAccumulator &errAccum);


void printDeviceArray(float *array_dev, int width, int height, int depth,
                      const char *name = NULL);


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
                    float *l1Error, float *l2Error, float *mse, float *pSNR) {

  WaveletCompressionParam param;
  param.binCount = binCount;
  param.quantAlg = quantAlg;
  param.thresholdValue = thresholdValue;
  param.waveletAlg = o->waveletAlg;
  param.isWaveletTransposeStandard = true;
  param.transformSteps = o->transformSteps;
  param.originalSize = o->originalData->size;

  *outputSizeBytes = 0;

  int count = o->count();

  float *inverseWaveletInput_dev;
  CUCHECK(cudaMalloc((void**)&inverseWaveletInput_dev, count * sizeof(float)));

  // If binCount is nonpositive, don't quantize and dequantize.
  // Just apply the threshold and reverse the wavelet transform.
  if (binCount <= 0) {

    // apply threshold
    // copy the data to a new array, zeroing out values less than threshold

    applyThresholdKernel<<<16,1024>>>
      (inverseWaveletInput_dev, (float*)o->transformedData->data_, count,
       thresholdValue);

  } else {  // binCount > 0

    // in one pass over the data:
    //   read one element
    //   quantize value
    //   increment appropriate histogram bin
    //   dequantize value
    //   write to new array

    /*
    quantHistDequantKernel<<<16,1024>>>
      (inverseWaveletInput_dev, (float*)transformedData.data_, count, 
    */ 

  }

  // compare to original
  ErrorAccumulator errAccum;

  float *tempData_dev;
  CUCHECK(cudaMalloc((void**)&tempData_dev, count * sizeof(float)));

  computeErrorRatesAfterDequantGPU
    (inverseWaveletInput_dev, o->transformedData->size,
     tempData_dev, param, (const unsigned char *)o->originalData->data_,
     errAccum);
     

  CUCHECK(cudaFree(tempData_dev));
  CUCHECK(cudaFree(inverseWaveletInput_dev));

  *l1Error = errAccum.getL1Error();
  *l2Error = errAccum.getL2Error();
  *mse = errAccum.getMeanSquaredError();
  *pSNR = errAccum.getPeakSignalToNoiseRatio();

  return true;
}
