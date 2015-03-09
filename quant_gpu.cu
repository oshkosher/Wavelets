#include "quant.h"
#include "quant_gpu.h"
#include "cudalloyds.h"
#include "test_compress_common.h"
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include "test_compress_gpu.h"

#define QUANTIZE_BLOCK_SIZE 1024
#define QUANTIZE_BLOCK_COUNT 64

#define MAX_CODEBOOK_SIZE 2048

// __constant__ float constCodebookArray[MAX_CODEBOOK_SIZE];

struct QuantLogFunctor {
  QuantLog uni;

  __host__ __device__ QuantLogFunctor(int binCount, float threshold,
                                      float maxVal) {
    uni.init(binCount, threshold, maxVal);
  }

  __host__ __device__ int operator() (float x) const {
    return uni.quant(x);
  }
};


struct DequantLogFunctor {
  QuantLog uni;

  __host__ __device__ DequantLogFunctor(int binCount, float threshold,
                                        float maxVal) {
    uni.init(binCount, threshold, maxVal);
  }

  __host__ __device__ float operator() (int x) const {
    return uni.dequant(x);
  }
};


struct QuantUniformFunctor {
  QuantUniform uni;

  __host__ __device__ QuantUniformFunctor(int binCount, float threshold,
                                          float maxVal) {
    uni.init(binCount, threshold, maxVal);
  }

  __host__ __device__ int operator() (float x) const {
    return uni.quant(x);
  }
};


struct DequantUniformFunctor {
  QuantUniform uni;

  __host__ __device__ DequantUniformFunctor(int binCount, float threshold,
                                            float maxVal) {
    uni.init(binCount, threshold, maxVal);
  }

  __host__ __device__ float operator() (int x) const {
    return uni.dequant(x);
  }
};


struct QuantCodebookFunctor {
  int boundaryCount;
  const float *boundaries_dev;

  static float *copyVectorToGPU(const vector<float> &v) {
    float *v_dev;
    size_t bytes = v.size()*sizeof(float);
    CUCHECK(cudaMalloc((void**)&v_dev, bytes));
    CUCHECK(cudaMemcpy(v_dev, v.data(), bytes, cudaMemcpyHostToDevice));
    return v_dev;
  }

  __host__ __device__ QuantCodebookFunctor(int boundaryCount_,
                                          const float *boundaries_dev_) {
    boundaryCount = boundaryCount_;
    boundaries_dev = boundaries_dev_;
  }

  __host__ __device__ int operator() (float x) const {
    return QuantCodebook::quant(x, boundaryCount, boundaries_dev);
  }
};


struct DequantCodebookFunctor {
  int binCount;
  const float *binValues;

  __host__ __device__ DequantCodebookFunctor(int binCount_,
                                             const float *binValues_) {
    binCount = binCount_;
    binValues = binValues_;
  }

  __host__ __device__ float operator() (int x) const {
    unsigned ux = x;
    if (ux >= binCount) {
      return 0.0f;
    } else {
      return binValues[ux];
    }
  }
};


bool quantizeGPU(int *outputData, const float *inputData, int count,
                 WaveletCompressionParam &param,
                 const float *nonzeroData_dev, int nonzeroCount,
                 float maxAbsVal, float minValue, float maxValue,
                 CudaTimer &quantizeTimer, int *zeroBin) {
  thrust::device_ptr<const float> inputStart(inputData), 
    inputEnd(inputData+count);
  thrust::device_ptr<int> outputStart(outputData);

  switch (param.quantAlg) {
  case QUANT_ALG_UNIFORM:
    {
      param.maxValue = maxAbsVal;
      quantizeTimer.start();
      QuantUniformFunctor q(param.binCount, param.thresholdValue,
                            param.maxValue);
      thrust::transform(inputStart, inputEnd, outputStart, q);
      quantizeTimer.end();
      break;
    }

  case QUANT_ALG_LOG:
    {
      param.maxValue = maxAbsVal;
      quantizeTimer.start();
      QuantLogFunctor q(param.binCount, param.thresholdValue, param.maxValue);
      thrust::transform(inputStart, inputEnd, outputStart, q);
      quantizeTimer.end();
      break;
    }

  case QUANT_ALG_LLOYD:
    {
      computeLloydQuantizationGPU(nonzeroData_dev, nonzeroCount, param.binCount,
                                  minValue, maxValue, param.thresholdValue,
                                  param.binBoundaries, param.binValues);
      quantizeTimer.start();
      float *boundaries_dev =
        QuantCodebookFunctor::copyVectorToGPU(param.binBoundaries);
      QuantCodebookFunctor q(param.binBoundaries.size(), boundaries_dev);
      thrust::transform(inputStart, inputEnd, outputStart, q);
      CUCHECK(cudaFree(boundaries_dev));
      quantizeTimer.end();
      break;
    }

  default:
    fprintf(stderr, "Quantization algorithm %d not found.\n",
            (int)param.quantAlg);
    return false;
  }

  // compute the bin number to which zero values map
  if (zeroBin) {
    Quantizer *quantizer = createQuantizer(param);
    *zeroBin = quantizer->quant(0);
    delete quantizer;
  }

  return true;
}


void computeLloydQuantizationGPU(const float *inputData, int count,
                                 int binCount, float minVal, float maxVal,
                                 float thresholdValue,
                                 std::vector<float> &quantBinBoundaries,
                                 std::vector<float> &quantBinValues) {
  
  // Make 'codebookSize' entries on either size of 0
  int codebookSize = (binCount-1) / 2;

  quantBinBoundaries.clear();
  quantBinValues.clear();

  // inputData is sorted, use it to get minVal and maxVal
  // Skip the last entry in inputData[] because it is often much larger than
  // the rest and skews the results
  float maxAbsVal;
  CUCHECK(cudaMemcpy(&maxAbsVal, inputData+count-2, sizeof(float),
                     cudaMemcpyDeviceToHost));
  assert(maxAbsVal > 0);

  vector<float> codebook;
  initialLloydCodebook(codebook, codebookSize, thresholdValue, maxAbsVal);

  /*
  printf("Before Lloyd\n");
  for (int i=0; i < codebookSize; i++) printf("%d) %f\n", i, codebook[i]);
  */

  // fine-tune the codebook and bin boundaries using Lloyd's algorithm.
  // This also applies the quantization to each value, writing the values
  // to quantizedData[].
  int lloydIters = 0;
  CudaTimer lloydTimer;
  lloydTimer.start();
  cudaLloyd(inputData, count-1, codebook.data(), (int)codebook.size(),
            DEFAULT_LLOYD_STOP_CRITERIA, true, &lloydIters);
  lloydTimer.end();
  lloydTimer.sync();
  if (!QUIET)
    printf("GPU lloyd %d iterations, %.3f ms\n", lloydIters, lloydTimer.time());
  
  /*  
  printf("After Lloyd\n");
  for (int i=0; i < codebookSize; i++) printf("%d) %f\n", i, codebook[i]);
  */

  setBinsFromCodebook(quantBinValues, quantBinBoundaries, binCount,
                      codebook, thresholdValue, minVal, maxVal);
}


// change data_dev from int[] to float[] in place
bool dequantizeGPU(float *outputData, const int *inputData,
                   int count, const WaveletCompressionParam &param) {

  CudaTimer timer("Dequantize");
  thrust::device_ptr<const int> inputStart(inputData), 
    inputEnd(inputData+count);
  thrust::device_ptr<float> outputStart(outputData);


  if (inputData == NULL)
    inputData = (const int*) outputData;
  
  timer.start();
  switch (param.quantAlg) {
  case QUANT_ALG_UNIFORM:
    {
      DequantUniformFunctor q(param.binCount, param.thresholdValue,
                              param.maxValue);
      thrust::transform(inputStart, inputEnd, outputStart, q);
      break;
    }
    
  case QUANT_ALG_LOG:
    {
      DequantLogFunctor q(param.binCount, param.thresholdValue, param.maxValue);
      thrust::transform(inputStart, inputEnd, outputStart, q);
      break;
    }

  case QUANT_ALG_LLOYD:
    {
      float *values_dev =
        QuantCodebookFunctor::copyVectorToGPU(param.binValues);
      DequantCodebookFunctor q(param.binCount, values_dev);
      thrust::transform(inputStart, inputEnd, outputStart, q);
      break;
    }

  default:
    fprintf(stderr, "Quantization algorithm %d not found.\n",
            (int)param.quantAlg);
    return false;
  }
  timer.end();
  if (!QUIET) {
    timer.sync();
    timer.print();
  }

  // copy the data to the CPU and print it all
  /*
  int *input = new int[count];
  float *output = new float[count];
  CUCHECK(cudaMemcpy(input, inputData, count*sizeof(int),
                     cudaMemcpyDeviceToHost));
  CUCHECK(cudaMemcpy(output, outputData, count*sizeof(float),
                     cudaMemcpyDeviceToHost));
  for (int i=0; i < count; i++)
    printf("%d) %d -> %f\n", i, input[i], output[i]);
  delete[] input;
  delete[] output;
  */
  return true;
}

