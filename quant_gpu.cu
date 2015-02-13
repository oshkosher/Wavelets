#include "quant.h"
#include "quant_gpu.h"
#include "cudalloyds.h"
#include "test_compress_common.h"


#define QUANTIZE_BLOCK_SIZE 1024
#define QUANTIZE_BLOCK_COUNT 64

#define MAX_CODEBOOK_SIZE 2048

// __constant__ float constCodebookArray[MAX_CODEBOOK_SIZE];

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


bool quantizeGPU(int *outputData, const float *inputData, int count,
                 WaveletCompressionParam &param,
                 const float *nonzeroData_dev, int nonzeroCount,
                 float maxAbsVal, float minValue, float maxValue,
                 CudaTimer &quantizeTimer) {
  switch (param.quantAlg) {
  case QUANT_ALG_UNIFORM:
    param.maxValue = maxAbsVal;
    quantizeTimer.start();
    quantUniformKernel<<<QUANTIZE_BLOCK_COUNT,QUANTIZE_BLOCK_SIZE>>>
      (outputData, inputData, count, param.binCount, param.thresholdValue,
       param.maxValue);
    quantizeTimer.end();
    break;

  case QUANT_ALG_LOG:
    param.maxValue = maxAbsVal;
    quantizeTimer.start();
    quantLogKernel<<<QUANTIZE_BLOCK_COUNT,QUANTIZE_BLOCK_SIZE>>>
      (outputData, inputData, count, param.binCount, param.thresholdValue,
       param.maxValue);
    quantizeTimer.end();
    break;

  case QUANT_ALG_LLOYD:
    computeLloydQuantizationGPU(nonzeroData_dev, nonzeroCount, param.binCount,
                                minValue, maxValue, param.thresholdValue,
                                param.binBoundaries, param.binValues);
    quantizeTimer.start();
    quantCodebookGPU(outputData, inputData, count, param.binBoundaries,
                     param.binCount);
    quantizeTimer.end();
    break;

  default:
    fprintf(stderr, "Quantization algorithm %d not found.\n",
            (int)param.quantAlg);
    return false;
  }

  return true;
}


template<class Q>
void __device__ quantKernel(const Q &quanter, int *output, const float *input,
                            int count) {

  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  while (idx < count) {
    float in = input[idx];
    output[idx] = quanter.quant(in);
    /*
    printf("[%5d of %d] %d,%d: %f -> %d\n",
           idx, count, blockIdx.x, threadIdx.x, in, output[idx]);
    */
    idx += blockDim.x * gridDim.x;
  }
}


void __global__ quantLogKernel(int *output, const float *input, int count,
                               int binCount, float threshold, float max) {
  __shared__ QuantLog quanter;

  if (threadIdx.x == 0) {
    quanter.init(binCount, threshold, max);
    // printf("GPU %.10g->%d\n", -0.009803920984, quanter.quant(-0.009803920984));
  }

  __syncthreads();

  quantKernel(quanter, output, input, count);
}


void __global__ quantUniformKernel(int *output, const float *input, int count,
                                   int binCount, float threshold, float max) {
  __shared__ QuantUniform quanter;

  if (threadIdx.x == 0) {
    quanter.init(binCount, threshold, max);
  }

  __syncthreads();

  quantKernel(quanter, output, input, count);
}

class QuantCodebookDev {
  int boundaryCount;
  const float *boundaries;
public:
  __device__ void init(int c, const float *b) {
    boundaryCount = c;
    boundaries = b;
  }

  __device__ int quant(float x) const {
    return QuantCodebook::quant(x, boundaryCount, boundaries);
  }
};

/*
class QuantCodebookConstDev {
  int boundaryCount;

public:
  __device__ void init(int c) {
    boundaryCount = c;
  }

  __device__ int quant(float x) const {
    return QuantCodebook::quant(x, boundaryCount, constCodebookArray);
  }
};
*/

void __global__ quantCodebookKernel(int *output, const float *input, int count, 
                                    const float *boundaries,
                                    int boundaryCount) {
  QuantCodebookDev quanter;
  quanter.init(boundaryCount, boundaries);
  quantKernel(quanter, output, input, count);
}


void quantCodebookGPU(int *output, const float *input, int dataCount,
                      const vector<float> &boundaries, int binCount) {

  // copy the boundaries array to the GPU
  size_t boundariesBytes = sizeof(float)*(binCount-1);
  float *boundaries_dev = NULL;
  CUCHECK(cudaMalloc((void**)&boundaries_dev, boundariesBytes));
  CUCHECK(cudaMemcpy(boundaries_dev, boundaries.data(), boundariesBytes,
                     cudaMemcpyHostToDevice));

  /*
  if (binCount > MAX_CODEBOOK_SIZE) {
    CUCHECK(cudaMalloc((void**)&boundaries_dev, boundariesBytes));
    CUCHECK(cudaMemcpy(boundaries_dev, boundaries.data(), boundariesBytes,
                       cudaMemcpyHostToDevice));
  } else {
    // this actually turned out slower
    CUCHECK(cudaMemcpyToSymbol(constCodebookArray, boundaries.data(),
                               boundariesBytes, 0, cudaMemcpyHostToDevice));
  }
  */

  quantCodebookKernel<<<QUANTIZE_BLOCK_COUNT,QUANTIZE_BLOCK_SIZE>>>
    (output, input, dataCount, boundaries_dev, binCount-1);

  if (boundaries_dev)
    CUCHECK(cudaFree(boundaries_dev));
}


template<class Q>
void __device__ dequantKernel(const Q &quanter, float *output,
                              const int *input, int count) {

  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  while (idx < count) {
    int in = input[idx];
    output[idx] = quanter.dequant(in);
    // printf("dequant %6d: %d %.2g\n", idx, in, output[idx]);
    idx += blockDim.x * gridDim.x;
  }
}


void __global__ dequantUniformKernel(float *output, const int *input, int count,
                                     int binCount, float threshold, float max) {
  __shared__ QuantUniform quanter;

  if (threadIdx.x == 0) {
    quanter.init(binCount, threshold, max);
  }

  __syncthreads();

  dequantKernel(quanter, output, input, count);
}


void __global__ dequantLogKernel(float *output, const int *input, int count,
                                 int binCount, float threshold, float max) {
  __shared__ QuantLog quanter;
  /*
  QuantLog quanter;
  quanter.init(binCount, threshold, max);
  */

  if (threadIdx.x == 0) {
    quanter.init(binCount, threshold, max);
    // printf("GPU %.10g->%d\n", -0.009803920984, quanter.quant(-0.009803920984));
  }

  __syncthreads();


  dequantKernel(quanter, output, input, count);
}


class DequantCodebookDev {
  const float *binValues;
  int binCount;
public:
  __device__ void init(const float *binValues_, int binCount_) {
    binValues = binValues_;
    binCount = binCount_;
  }

  __device__ float dequant(int i) const {
    assert(i >= 0 && i < binCount);
    return binValues[i];
  }
};


void __global__ dequantCodebookKernel(float *output, const int *input,
                                      int count,
                                      const float *binValues, int binCount) {
  DequantCodebookDev dequanter;
  dequanter.init(binValues, binCount);

  dequantKernel(dequanter, output, input, count);
}


void dequantCodebookGPU(float *result_dev, const int *input_dev, int count,
                        const vector<float> &binValues) {
  float *binValues_dev;
  size_t binValuesBytes = binValues.size() * sizeof(float);
  CUCHECK(cudaMalloc((void**)&binValues_dev, binValuesBytes));
  CUCHECK(cudaMemcpy(binValues_dev, binValues.data(), binValuesBytes,
                     cudaMemcpyHostToDevice));

  dequantCodebookKernel<<<QUANTIZE_BLOCK_COUNT,QUANTIZE_BLOCK_SIZE>>>
    (result_dev, input_dev, count, binValues_dev, (int)binValues.size());

  CUCHECK(cudaFree(binValues_dev));
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
  CUCHECK(cudaThreadSynchronize());
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
bool dequantizeGPU(float *result_dev, const int *input_dev,
                   int count, const WaveletCompressionParam &param) {

  CudaTimer timer("Dequantize");

  if (input_dev == NULL)
    input_dev = (const int*) result_dev;
  
  timer.start();
  switch (param.quantAlg) {
  case QUANT_ALG_UNIFORM:
    dequantUniformKernel<<<QUANTIZE_BLOCK_COUNT,QUANTIZE_BLOCK_SIZE>>>
      (result_dev, input_dev, count,
       param.binCount, param.thresholdValue, param.maxValue);
    break;

  case QUANT_ALG_LOG:
    dequantLogKernel<<<QUANTIZE_BLOCK_COUNT,QUANTIZE_BLOCK_SIZE>>>
      (result_dev, input_dev, count,
       param.binCount, param.thresholdValue, param.maxValue);
    break;

  case QUANT_ALG_LLOYD:
    dequantCodebookGPU(result_dev, input_dev, count, param.binValues);
    break;

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
  CUCHECK(cudaMemcpy(input, input_dev, count*sizeof(int),
                     cudaMemcpyDeviceToHost));
  CUCHECK(cudaMemcpy(output, result_dev, count*sizeof(float),
                     cudaMemcpyDeviceToHost));
  for (int i=0; i < count; i++)
    printf("%d) %d -> %f\n", i, input[i], output[i]);
  delete[] input;
  delete[] output;
  */
  return true;
}

