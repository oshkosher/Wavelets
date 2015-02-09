#include "cuda.h"
#include <thrust/execution_policy.h>
#include <thrust/functional.h>
#include <thrust/transform.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/binary_search.h>
#include <thrust/sort.h>
#include "test_compress_common.h"
#include "dquant_log_cpu.h"
#include "dquant_unif_cpu.h"
#include "dwt_cpu.h"
#include "dwt_gpu.h"
#include "nixtimer.h"
#include "quant_count.h"
#include "quant_log_cpu.h"
#include "quant_unif_cpu.h"
#include "thresh_cpu.h"
#include "cucheck.h"
#include "cuda_timer.h"
#include "quant.h"
#include "cudalloyds.h"
#include "CUDA/wavelet/wavelet.h"

#define QUANTIZE_BLOCK_SIZE 1024
#define QUANTIZE_BLOCK_COUNT 64

// #define DO_DETAILED_CHECKS

bool compressFile(const char *inputFile, const char *outputFile, Options &opt);
bool decompressFile(const char *inputFile, const char *outputFile,
		    Options &opt);

bool check1WaveletTransform(const float *data_dev, const CubeFloat *dataBefore,
                            const WaveletCompressionParam &param,
                            bool verbose);

bool check2Sorted(const float *data_dev, const float *orig_data_dev,
                  int dataCount);

bool check3Quantized(const int *quantizedData_dev, int dataCount,
                     const float *unquantizedData_dev,
                     const WaveletCompressionParam &param,
                     scu_wavelet::int3 size,
                     float thresholdValue, float max, bool verbose);

bool waveletTransformGPU(float *data_d, float *tmpData_d,
                         scu_wavelet::int3 &size, 
                         const WaveletCompressionParam &param,
                         bool isInverse, CudaTimer *waveletTimer,
                         CudaTimer *transposeTimer);

void sortAbsAndFindLimitsGPU(float *dest_d, const float *src_d, int count,
                             float *minValue, float *maxValue);

void computeLloydQuantizationGPU(const float *inputData, int count,
                                 int binCount, float minVal, float maxVal,
                                 float thresholdValue,
                                 std::vector<float> &quantBinBoundaries,
                                 std::vector<float> &quantBinValues);

bool computeErrorRatesGPU(int *quantizedData_dev,
                          scu_wavelet::int3 size,
                          float *tempData_dev,
                          const WaveletCompressionParam &param,
                          const unsigned char *inputData_dev,
                          float *meanSqErr, float *peakSNR);

void computeErrorRatesAfterDequantGPU
(float *data_dev, scu_wavelet::int3 size, float *tempData_dev,
 const WaveletCompressionParam &param,
 const unsigned char *inputData_dev, ErrorAccumulator &errAccum);


// If input_dev is NULL, use (int*)result_dev as the input.
bool dequantizeGPU(float *result_dev, const int *input_dev,
                   int count, const WaveletCompressionParam &param);

void __global__ quantUniformKernel(void *data, int count,
                                   int binCount, float threshold, float max);
void __global__ quantLogKernel(void *data, int count,
                               int binCount, float threshold, float max);
void __global__ quantCodebookKernel(void *data, int count, 
                                    const float *boundaries, int boundaryCount);
void quantCodebookGPU(void *data_dev, int dataCount,
                      const vector<float> &boundaries, int binCount);

void __global__ dequantUniformKernel(float *output, const int *input, int count,
                                     int binCount, float threshold, float max);
void __global__ dequantLogKernel(float *data, const int *input, int count,
                                 int binCount, float threshold, float max);
void __global__ dequantCodebookKernel(float *output, const int *input,
                                      int count,
                                      const float *binValues, int binCount);
void dequantCodebookGPU(float *result_dev, const int *input_dev, int count,
                        const vector<float> &binValues);

inline bool isClose(float a, float b) {
  float diff = fabsf(a-b);
  return ((a==0 || b==0) ? diff : diff/a) < 0.001;
}

inline bool isClose(int a, int b) {
  return abs(a-b) <= 1;
}

// return the number of mismatches
int compareArrays(const float *a, const float *b, int count);
int compareArrays(const int *a, const int *b, int count);

// debug output
void printArray(float *array, int width, int height, int depth,
                const char *name = NULL);
void printArray(int *array, int width, int height, int depth,
                const char *name = NULL);
void printDeviceArray(float *array_dev, int width, int height, int depth,
                      const char *name = NULL);
void printDeviceArray(float *array_dev, scu_wavelet::int3 size,
                      const char *name = NULL);

// struct AbsFunctor : public thrust::unary_function<float,float> {
struct AbsFunctor {
  __host__ __device__ float operator() (float x) const {
    return fabsf(x);
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
  

/**
   Do a Haar wavelet transform on data_dev[]. The array is already on
   the GPU. The data is a square size*size and size is a power of 2.
   data_dev[] contains the input data and will be overwritten with the
   output data. 

   tmp_data_dev[] - used as temporary storage during the computation.
     Must be at least as large as data_dev[].
   stepCount - # of wavelet steps
   isInverse - iff true, do an inverse transformation
   isStandardTranspose -
      if true: do all xform steps for one direction
               transpose
               do all xform steps for other direction
               transpose

      if false: do a single xform step for one direction
                transpose
                do a single xform step for other direction
                transpose
                repeat for each step
*/
bool wavelet_transform_gpu(float *data_dev, float *tmp_data_dev,
                           int size, int stepCount, bool isInverse,
                           bool isStandardTranspose, cudaStream_t stream=0);


int main(int argc, char **argv) {

  Options opt;
  int nextArg;

  if (!parseOptions(argc, argv, opt, nextArg)) return 1;
  
  // set global variable to enable/disable status output
  QUIET = opt.quiet;

  if (argc - nextArg != 2) printHelp();

  const char *inputFile = argv[nextArg++];
  const char *outputFile = argv[nextArg++];
  bool result;

  if (opt.doCompress) {
    result = compressFile(inputFile, outputFile, opt);
  } else {
    result = decompressFile(inputFile, outputFile, opt);
  }

  // deallocate static protobuf data
  google::protobuf::ShutdownProtobufLibrary();

  if (result == false) return 1;

  return 0;
}

bool compressFile(const char *inputFile, const char *outputFile,
                  Options &opt) {

  WaveletCompressionParam param = opt.param;
  Cube inputData;
  CubeFloat data;
  CubeInt quantizedData;
  CubeletStreamReader cubeletStream;
  double firstStartTime, startTime;
  CudaTimer copyToGPUTimer("Copy to GPU"), waveletTimer("Wavelet transform"),
    transposeTimer("Transpose"),
    /* sortTimer("Sort"), */ quantizeTimer("Quantize"), 
    copyFromGPUTimer("Copy from GPU");
  // FILE *f = fopen("marker","w"); fprintf(f, "marker\n"); fclose(f);
  firstStartTime = NixTimer::time();

  // on CPU: read the data into memory, convert to float, pad
  // optional: convert to float on the GPU
  // optional: pad on the GPU
  if (!readData(inputFile, &inputData, &data, &cubeletStream, &param))
    return false;

  // if not there already, copy to GPU
  // XXX use a stream for asynchronous operation

  // allocate memory on the GPU for two copies of the data
  unsigned char *inputData_dev = NULL;
  float *data1_dev, *data2_dev;
  size_t dataCount = data.count();
  size_t dataCountBytes = dataCount * sizeof(float);
  startTime = NixTimer::time();
  CUCHECK(cudaMalloc((void**)&data1_dev, dataCountBytes));
  CUCHECK(cudaMalloc((void**)&data2_dev, dataCountBytes));

  // if the input data is bytes, we can compute its error rate later.
  // copy it to the GPU so we can do that
  if (inputData.datatype == WAVELET_DATA_UINT8) {
    CUCHECK(cudaMalloc((void**)&inputData_dev, dataCount));
  }

  if (!QUIET) {
    printf("CUDA initialize: %.2f ms\n", 1000*(NixTimer::time() - startTime));
    fflush(stdout);
  }

  // copy to GPU
  copyToGPUTimer.start();
  CUCHECK(cudaMemcpy(data1_dev, data.data(), dataCountBytes,
                     cudaMemcpyHostToDevice));
  if (inputData_dev)
    CUCHECK(cudaMemcpy(inputData_dev, inputData.data_, dataCount,
                       cudaMemcpyHostToDevice));
  copyToGPUTimer.end();

  // size of the data as it exists on the device
  scu_wavelet::int3 deviceSize = data.size;

  // Wavelet transform
  if (opt.verbose) printDeviceArray(data1_dev, deviceSize.x, deviceSize.y,
                                    deviceSize.z, "Before fwd");
  waveletTransformGPU(data1_dev, data2_dev, deviceSize, param, false,
                      &waveletTimer, &transposeTimer);
  if (opt.verbose) printDeviceArray(data1_dev, deviceSize.x, deviceSize.y,
                                    deviceSize.z, "After fwd");

  if (!QUIET) {
    // sync to get the timer values
    CUCHECK(cudaThreadSynchronize());
    copyToGPUTimer.print();
    printf("Wavelet transform: %.3f ms (transpose %.3f ms)\n",
           waveletTimer.time() + transposeTimer.time(),
           transposeTimer.time());
  }

  // just for testing--copy the data back to the CPU to see if it's correct
#ifdef DO_DETAILED_CHECKS
  check1WaveletTransform(data1_dev, &data, param, opt.verbose);
#endif

  // sort the absolute values and compute the positive and negative values
  // with the greatest magnitude
  // input is data1_dev, sorted output is data2_dev
  float minValue, maxValue;
  sortAbsAndFindLimitsGPU(data2_dev, data1_dev, dataCount,
                          &minValue, &maxValue);

  // copy the data to data2_dev to be sorted
  // CUCHECK(cudaMemcpy(data2_dev, data1_dev, dataCountBytes,
  // cudaMemcpyDeviceToDevice));

  // sort the absolute values on the GPU
  thrust::device_ptr<float>
    data1_start(data1_dev),
    data1_end(data1_dev + dataCount),
    data2_start(data2_dev),
    data2_end(data2_dev + dataCount);


  /* thrust::sort works with pointers on the device (like data2_dev), 
     as long as you give it "thrust::device" as the first argument.
     thrust::transform does not. It requires thrust::device_ptr objects
     (like data1_start). thrust::sort will work with those too.
  */
  // XXX try out sorting with a transform_iterator

  // make a copy of the data, applying abs() to each value, into data2,
  // then sort that copy

  /*
  AbsFunctor absFunctor;
  sortTimer.start();
  thrust::transform(data1_start, data1_end, data2_start, absFunctor);
  thrust::sort(data2_start, data2_start + dataCount);
  sortTimer.end();
  */

#ifdef DO_DETAILED_CHECKS
  CUCHECK(cudaThreadSynchronize());
  check2Sorted(data2_dev, data1_dev, dataCount);
#endif

  // data1_dev now contains the data, and data2_dev contains the sorted
  // absolute values of the data

  // fetch percentage threshold values from GPU

  // send back the min, minimum nonzero, and maximum

  int thresholdOffset;
  float max = std::max(fabsf(minValue), fabsf(maxValue));

  // if you omit thrust::less<float>(), it defaults to an integer comparison,
  // and looks for the first value greater than or equal to 1
  /*
  float minNonzero;
  size_t minNonzeroOffset = 
    thrust::upper_bound(data2_start, data2_end, 0, thrust::less<float>())
    - data2_start;
  CUCHECK(cudaMemcpy(&minNonzero, data2_dev + minNonzeroOffset,
                     sizeof(float), cudaMemcpyDeviceToHost));
  */

  // get the threshold value
  if (param.thresholdFraction <= 0) {
    param.thresholdValue = MIN_THRESHOLD_VALUE;
    thresholdOffset = 0;
  } else if (param.thresholdFraction >= 1) {
    param.thresholdValue = max;
    thresholdOffset = dataCount;
  } else {
    thresholdOffset = (int) (param.thresholdFraction * dataCount);
    // printf("threshold offest = %d\n", thresholdOffset);
    CUCHECK(cudaMemcpy(&param.thresholdValue, data2_dev + thresholdOffset,
                       sizeof(float), cudaMemcpyDeviceToHost));
    if (param.thresholdValue == 0)
      param.thresholdValue = MIN_THRESHOLD_VALUE;
  }

  if (!QUIET) {
    printf("min = %f, max = %f, threshold = %.10g\n",
           minValue, maxValue, param.thresholdValue);
    fflush(stdout);
  }

  // choose threshold and bin count on CPU

  // quantize on GPU
#ifdef DO_DETAILED_CHECKS
  float *unquantizedCopy = NULL;
  CUCHECK(cudaMalloc((void**)&unquantizedCopy, dataCountBytes));
  CUCHECK(cudaMemcpy(unquantizedCopy, data1_dev, dataCountBytes,
                     cudaMemcpyDeviceToDevice));
#endif

  // overwrite the data in data1_dev with quantized integers

  quantizeTimer.start();
  switch (param.quantAlg) {
  case QUANT_ALG_UNIFORM:
    param.maxValue = max;
    quantUniformKernel<<<QUANTIZE_BLOCK_COUNT,QUANTIZE_BLOCK_SIZE>>>
      (data1_dev, dataCount,
       param.binCount, param.thresholdValue, param.maxValue);
    break;

  case QUANT_ALG_LOG:
    param.maxValue = max;
    quantLogKernel<<<QUANTIZE_BLOCK_COUNT,QUANTIZE_BLOCK_SIZE>>>
      (data1_dev, dataCount,
       param.binCount, param.thresholdValue, param.maxValue);
    break;

  case QUANT_ALG_LLOYD:
    const float *nonzeroData_d;
    int nonzeroCount;
    nonzeroData_d = data2_dev + thresholdOffset;
    nonzeroCount = dataCount - thresholdOffset;
    computeLloydQuantizationGPU(nonzeroData_d, nonzeroCount, param.binCount,
                                minValue, maxValue, param.thresholdValue,
                                param.binBoundaries, param.binValues);
    quantCodebookGPU(data1_dev, dataCount, param.binBoundaries, param.binCount);
    break;

  default:
    fprintf(stderr, "Quantization algorithm %d not found.\n",
            (int)param.quantAlg);
    return false;
  }
  quantizeTimer.end();

#ifdef DO_DETAILED_CHECKS
  if (unquantizedCopy)
    check3Quantized((int*)data1_dev, dataCount, unquantizedCopy,
                    param, deviceSize, param.thresholdValue, max, opt.verbose);
#endif

  // this is needed to get the final timer value
  CUCHECK(cudaThreadSynchronize());

  if (!QUIET) {
    quantizeTimer.print();
    fflush(stdout);
  }

  // compute error rates
  if (opt.doComputeError) {
    float meanSqErr, peakSNR;
    // quantizedData.print("quantized before computing error rates");
    if (computeErrorRatesGPU((int*)data1_dev, deviceSize, data2_dev,
                             param, inputData_dev, &meanSqErr, &peakSNR)) {
      printf("Mean squared error: %.3f, peak SNR: %.3f\n", meanSqErr, peakSNR);
      fflush(stdout);
    }
  }

  // copy the data back to the CPU
  copyFromGPUTimer.start();
  quantizedData.size = deviceSize;
  quantizedData.allocate();
  CUCHECK(cudaMemcpy(quantizedData.data(), data1_dev, dataCountBytes,
                     cudaMemcpyDeviceToHost));
  copyFromGPUTimer.end();
  if (!QUIET) {
    CUCHECK(cudaThreadSynchronize());
    copyFromGPUTimer.print();
  }

  // Huffman code and output

  // write the quantized data to a file
  // for now, make a new file just for this cubelet
  startTime = NixTimer::time();
  quantizedData.param = param;
  CubeletStreamWriter cubeletWriter;
  if (!cubeletWriter.open(outputFile)) return false;
  if (!writeQuantData(cubeletWriter, &quantizedData, opt))
    return false;
  cubeletWriter.close();
  
  if (!QUIET) {
    printf("Write data file: %.2f ms\n", (NixTimer::time() - startTime)*1000);
    printf("Total: %.2f ms\n", (NixTimer::time() - firstStartTime)*1000);
  }

  CUCHECK(cudaFree(data1_dev));
  CUCHECK(cudaFree(data2_dev));
  if (inputData_dev)
    CUCHECK(cudaFree(inputData_dev));

  return true;
}

bool decompressFile(const char *inputFile, const char *outputFile,
		    Options &opt) {

  /*
  FileData f;
  float *data;
  int width, height;
  double firstStartTime, startTime, elapsed;

  firstStartTime = startTime = NixTimer::time();

  if (!readQuantDataProtoBuf(inputFile, f)) return false;

  elapsed = NixTimer::time() - startTime;
  printf("Read %dx%d data file: %.2f ms\n", f.width, f.height,
	 (NixTimer::time() - startTime) * 1000);

  data = f.data;
  width = f.width;
  height = f.height;

  assert(data != NULL && width > 0 && height > 0);
  assert(width == height);
  int size = width;

  // de-quantize the data
  startTime = NixTimer::time();
  switch (f.quantizeAlgorithm) {
  case QUANT_ALG_UNIFORM:
    dquant_unif_cpu(size, data, f.quantizeBits, f.threshold, f.quantMaxVal);
    break;
  case QUANT_ALG_LOG:
    dquant_log_cpu(size, data, f.quantizeBits, f.threshold, f.quantMaxVal);
    break;
  case QUANT_ALG_COUNT:
    // check the size of the codebook
    dequant_codebook_array(f.quantBinValues, size*size, data);
    break;
  case QUANT_ALG_LLOYD:
    fprintf(stderr, "Lloyd's algorithm not integrated yet.\n");
    return false;
  default:
    fprintf(stderr, "Quantization algorithm %d not found.\n",
            (int)f.quantizeAlgorithm);
    return false;
  }

  printf("Dequantize: %.2f ms\n", (NixTimer::time() - startTime) * 1000);

  // perform inverse wavelet transform
  elapsed = haar_2d(size, data, true, f.waveletSteps);
  printf("Wavelet inverse transform: %.2f ms\n", elapsed);

  // write the reconstructed data
  startTime = NixTimer::time();
  if (!writeDataFile(outputFile, data, size, size, true)) return false;
  printf("Write file: %.2f ms\n", (NixTimer::time() - startTime) * 1000);

  printf("Total: %.2f ms\n", (NixTimer::time() - firstStartTime) * 1000);
  */
  return true;
}


bool waveletTransformGPU(float *data_d, float *tmpData_d,
                         scu_wavelet::int3 &size, 
                         const WaveletCompressionParam &param,
                         bool isInverse, CudaTimer *waveletTimer,
                         CudaTimer *transposeTimer) {
  switch (param.waveletAlg) {
  case WAVELET_HAAR:
    haar_3d_cuda(data_d, tmpData_d, size, param.transformSteps,
                 isInverse, waveletTimer, transposeTimer);
    return true;

  case WAVELET_CDF97:
    cdf97_3d_cuda(data_d, tmpData_d, size, param.transformSteps,
                  isInverse, waveletTimer, transposeTimer);
    return true;
    /*
    {
      float filter[9] = {
        .037828455506995,
        -.02384946501938,
        -.110624404418420,
        .377402855612650,
        .85269867900940,
        .377402855612650,
        -.110624404418420,
        -.02384946501938,
        .037828455506995
      };
      if (waveletTimer) waveletTimer->start();
      setUpFilter(filter);
      wavelet_cuda_3d_fwd(data_d, size.x, size.y, size.z,
                          param.transformSteps.x, param.transformSteps.y,
                          param.transformSteps.z, true);
      if (waveletTimer) waveletTimer->end();
    }
    return true;
    */


  default:
    fprintf(stderr, "Unrecognized wavelet id %d.\n", param.waveletAlg);
    return false;
  }
}


__global__ void copyAbsoluteAndMinMaxKernel(float *dest, const float *src,
                                            int count, float minMax[2]) {

  __shared__ float minShared, maxShared;
  float min = src[0];
  float max = min;

  if (threadIdx.x == 0) {
    minShared = maxShared = min;
  }
  
  __syncthreads();

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  while (idx < count) {
    float in = src[idx];
    if (in < min) min = in;
    if (in > max) max = in;
    
    dest[idx] = fabsf(in);

    idx += blockDim.x * gridDim.x;
  }

  WaveletAtomic::min(&minShared, min);
  WaveletAtomic::max(&maxShared, max);

  __syncthreads();

  if (threadIdx.x == 0) {
    WaveletAtomic::min(&minMax[0], minShared);
    WaveletAtomic::max(&minMax[1], maxShared);
  }
}
  

void sortAbsAndFindLimitsGPU(float *dest_d, const float *src_d, int count,
                             float *minValue, float *maxValue) {
  
  float *minMax_d;
  CUCHECK(cudaMalloc((void**)&minMax_d, sizeof(float)*2));
  CUCHECK(cudaMemcpy(minMax_d, src_d, sizeof(float),
                     cudaMemcpyDeviceToDevice));
  CUCHECK(cudaMemcpy(minMax_d+1, src_d, sizeof(float),
                     cudaMemcpyDeviceToDevice));

  CudaTimer absTimer("Absolute value"), sortTimer("Sort");
  absTimer.start();
  copyAbsoluteAndMinMaxKernel<<<4,1024>>>(dest_d, src_d, count, minMax_d);
  absTimer.end();

  sortTimer.start();
  thrust::sort(thrust::device, dest_d, dest_d + count);
  sortTimer.end();

  float minMax[2];
  CUCHECK(cudaMemcpy(minMax, minMax_d, sizeof(float)*2,
                     cudaMemcpyDeviceToHost));
  *minValue = minMax[0];
  *maxValue = minMax[1];
  CUCHECK(cudaFree(minMax_d));

  if (!QUIET) {
    absTimer.print();
    sortTimer.print();
  }
}


template<class Q>
void __device__ quantKernel(const Q &quanter, void *data, int count) {
  const float *dataFloat = (const float*) data;
  int *dataInt = (int*) data;

  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  while (idx < count) {
    float in = dataFloat[idx];
    dataInt[idx] = quanter.quant(in);
    // printf("[%6d] (%d,%d) %f -> %d\n", idx, blockIdx.x, threadIdx.x, in, dataInt[idx]);
    idx += blockDim.x * gridDim.x;
  }
}


void __global__ quantLogKernel(void *data, int count,
                               int binCount, float threshold, float max) {
  __shared__ QuantLog quanter;

  if (threadIdx.x == 0) {
    quanter.init(binCount, threshold, max);
    // printf("GPU %.10g->%d\n", -0.009803920984, quanter.quant(-0.009803920984));
  }

  __syncthreads();

  quantKernel(quanter, data, count);
}


void __global__ quantUniformKernel(void *data, int count,
                                   int binCount, float threshold, float max) {
  __shared__ QuantUniform quanter;

  if (threadIdx.x == 0) {
    quanter.init(binCount, threshold, max);
  }

  __syncthreads();

  quantKernel(quanter, data, count);
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

void __global__ quantCodebookKernel(void *data, int count, 
                                    const float *boundaries,
                                    int boundaryCount) {
  QuantCodebookDev quanter;
  quanter.init(boundaryCount, boundaries);
  
  quantKernel(quanter, data, count);
}


void quantCodebookGPU(void *data_dev, int dataCount,
                      const vector<float> &boundaries, int binCount) {

  // copy the boundaries array to the GPU
  float *boundaries_dev;
  size_t boundariesBytes = sizeof(float)*(binCount-1);
  CUCHECK(cudaMalloc((void**)&boundaries_dev, boundariesBytes));
  CUCHECK(cudaMemcpy(boundaries_dev, boundaries.data(), boundariesBytes,
                     cudaMemcpyHostToDevice));

  quantCodebookKernel<<<QUANTIZE_BLOCK_COUNT,QUANTIZE_BLOCK_SIZE>>>
    (data_dev, dataCount, boundaries_dev, binCount-1);

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
  


// return the number of mismatches
int compareArrays(const float *a, const float *b, int count) {
  int mismatches = 0;
  for (int i=0; i < count; i++) {
    // printf("[%6d] %g %g\n", i, a[i], b[i]);
    if (!isClose(a[i], b[i])) {
      mismatches++;
      // printf("mismatch at %d: %g  %g\n", i, a[i], b[i]);
    }
  }
  return mismatches;
}

// return the number of mismatches
int compareArrays(const int *a, const int *b, int count) {
  int mismatches = 0;
  for (int i=0; i < count; i++) {
    if (!isClose(a[i], b[i])) {
      mismatches++;
      // printf("mismatch at %d: %d  %d\n", i, a[i], b[i]);
    }
  }
  return mismatches;
}


void printArray(float *array, int width, int height, int depth, 
                const char *name) {
  if (name) printf("%s\n", name);
  for (int level=0; level < depth; level++) {
    printf("z=%d\n", level);
    for (int row=0; row < height; row++) {
      for (int col=0; col < width; col++) {
        printf("%8.4f ", array[(level*height + row)*width + col]);
      }
      putchar('\n');
    }
    putchar('\n');
  }
  putchar('\n');
}


void printArray(int *array, int width, int height, int depth, 
                const char *name) {
  if (name) printf("%s\n", name);
  for (int level=0; level < depth; level++) {
    printf("z=%d\n", level);
    for (int row=0; row < height; row++) {
      for (int col=0; col < width; col++) {
        printf("%5d ", array[(level*height + row)*width + col]);
      }
      putchar('\n');
    }
    putchar('\n');
  }
  putchar('\n');
}


void printDeviceArray(float *array_dev, scu_wavelet::int3 size,
                      const char *name) {
  printDeviceArray(array_dev, size.x, size.y, size.z, name);
}

void printDeviceArray(float *array_dev, int width, int height, int depth, 
                      const char *name) {

  float *array = new float[width*height*depth];
  CUCHECK(cudaMemcpy(array, array_dev, sizeof(float)*width*height*depth,
                     cudaMemcpyDeviceToHost));

  printArray(array, width, height, depth, name);

  delete[] array;
}


bool check1WaveletTransform(const float *data_dev, const CubeFloat *dataBefore,
                            const WaveletCompressionParam &param,
                            bool verbose) {
  
  // transform the data on the CPU
  CubeFloat result;
  result.size = dataBefore->size;
  result.allocate();
  result.copyFrom(*dataBefore);
  int count = dataBefore->count();

  // dataBefore->print("Before");

  // printArray(result.data(), result.size.x, result.size.y, result.size.z, "Before");
             
  
  // writeDataFile("check1a.data", data, width, height, true);
  // haar_3d(&result, param.transformSteps, false,
  // param.isWaveletTransposeStandard);
  // writeDataFile("check1b.data", data, width, height, true);
  waveletTransform(result, param, false, verbose);

  // copy the data from the GPU to the CPU
  float *dataTransformedByGPU = new float[count];
  CUCHECK(cudaMemcpy(dataTransformedByGPU, data_dev, sizeof(float) * count,
                     cudaMemcpyDeviceToHost));
  // writeDataFile("check1c.data", dataTransformedByGPU, width, height, true);

  // printArray(result.data(), result.size.x, result.size.y, result.size.z, "After CPU");

  if (verbose)
    printArray(dataTransformedByGPU, result.size.x, result.size.y,
               result.size.z, "After GPU");

  // printf("CPU.\n");
  // dataBefore->print();
  /*
  for (int i=0; i < dataBefore->size.x; i++) {
    printf("%d)\t%.4g\n", i, result.get(i,0,0));
  }
  */

  // printArray(dataTransformedByGPU, result.size.x, result.size.y, result.size.z, "GPU.");
             
             
  /*
  for (int i=0; i < dataBefore->size.x; i++) {
    printf("%d)\t%.4g\n", i, dataTransformedByGPU[i]);
  }
  */

  /*
  printf("After.\n");
  for (int i=0; i < count; i++) {
    printf("%d)\t%g\t%g\n", i, result.get(i,0,0), dataTransformedByGPU[i]);
  }
  */

  // compare the data
  int mismatches = compareArrays(result.data(), dataTransformedByGPU,
                                 count);

  printf("Check 1: wavelet transform %d / %d mismatches\n",
         mismatches, count);

  if (mismatches > 0) {
    // printArray(result.data(), result.size.x, result.size.y, result.size.z, "CPU");
    // printArray(dataTransformedByGPU, result.size.x, result.size.y, result.size.z, "GPU");
  }

  delete[] dataTransformedByGPU;

  return mismatches == 0;
}


bool check2Sorted(const float *data_dev, const float *orig_data_dev,
                  int dataCount) {

  float *cpuSorted = new float[dataCount];
  CUCHECK(cudaMemcpy(cpuSorted, orig_data_dev, sizeof(float) * dataCount,
                     cudaMemcpyDeviceToHost));
  for (int i=0; i < dataCount; i++) cpuSorted[i] = fabsf(cpuSorted[i]);
  std::sort(cpuSorted, cpuSorted + dataCount);

  float *gpuSorted = new float[dataCount];
  CUCHECK(cudaMemcpy(gpuSorted, data_dev, sizeof(float) * dataCount,
                     cudaMemcpyDeviceToHost));
  
  int mismatches = compareArrays(cpuSorted, gpuSorted, dataCount);

  printf("Check 2: sort %d / %d mismatches\n", mismatches, dataCount);

  /* for (int i=0; i < dataCount; i++) {
    printf("%.10g\t%.10g\n", cpuSorted[i], gpuSorted[i]);
    } */

  delete[] cpuSorted;
  delete[] gpuSorted;

  return mismatches == 0;
}


bool check3Quantized(const int *quantizedData_dev, int count,
                     const float *unquantizedData_dev,
                     const WaveletCompressionParam &param,
                     scu_wavelet::int3 size,
                     float thresholdValue, float max, bool verbose) {

  float *unquantized = new float[count];
  CUCHECK(cudaMemcpy(unquantized, unquantizedData_dev, count * sizeof(float),
                     cudaMemcpyDeviceToHost));

  int *gpuQuant = new int[count];
  CUCHECK(cudaMemcpy(gpuQuant, quantizedData_dev, count * sizeof(int),
                     cudaMemcpyDeviceToHost));

  Quantizer *quanter = createQuantizer(param);
  int *cpuQuant = new int[count];
  quanter->quantizeRow(unquantized, cpuQuant, count);

  // print dense quantized data
  if (verbose) {
    printArray(cpuQuant, size.x, size.y, size.z, "CPU quantized");
    printArray(gpuQuant, size.x, size.y, size.z, "GPU quantized");
  }

  int mismatches = compareArrays(cpuQuant, gpuQuant, count);
  int notEqual = 0;
  for (int i=0; i < count; i++) {
    // if (verbose) printf("[%d] %.10g  %d %d\n", i, unquantized[i], cpuQuant[i], gpuQuant[i]);
                        
    if (cpuQuant[i] != gpuQuant[i]) {
      notEqual++;
      if (verbose)
        printf("[%d] %.10g cpu %d gpu %d\n", i, unquantized[i], cpuQuant[i], gpuQuant[i]);
    }
  }

  printf("Check 3: quant %d / %d not close, %d not equal\n", mismatches, count, notEqual);


  /*
  for (int i=0; i < count; i++) {
    printf("%d\t%f\n", quantizedData[i], unquantizedData[i]);
  }
  */
  delete[] gpuQuant;
  delete[] cpuQuant;

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


// Size is transformed.
bool computeErrorRatesGPU(int *quantizedData_dev, scu_wavelet::int3 size,
                          float *tempData_dev,
                          const WaveletCompressionParam &param,
                          const unsigned char *inputData_dev,
                          float *meanSqErr, float *peakSNR) {

  float *data_dev;
  int count = size.count(), byteCount = size.count() * sizeof(float);

  CUCHECK(cudaMalloc((void**)&data_dev, byteCount));

  // dequantize quantizedData_dev into data_dev
  if (!dequantizeGPU(data_dev, quantizedData_dev, count, param)) return false;

  // printDeviceArray(data_dev, size.x, size.y, size.z, "Dequantized");

  ErrorAccumulator errAccum;

  computeErrorRatesAfterDequantGPU(data_dev, size, tempData_dev, param,
                                   inputData_dev, errAccum);
  
  *meanSqErr = errAccum.getMeanSquaredError();
  *peakSNR = errAccum.getPeakSignalToNoiseRatio();


  CUCHECK(cudaFree(data_dev));

  return true;
}


// errSums[0] = sumDiff
// errSums[1] = sumDiffSquared
__global__ void computeErrorKernel(float errSums[2],
                                   const float *data, 
                                   const unsigned char *orig,
                                   scu_wavelet::int3 dataSize,
                                   scu_wavelet::int3 origSize) {

  unsigned long long sumDiff = 0, sumDiffSquared = 0;
  __shared__ unsigned long long sumDiffShared, sumDiffSquaredShared;

  if (threadIdx.x == 0 && threadIdx.y == 0) {
    sumDiffShared = sumDiffSquaredShared = 0;
  }
 
  // sync so all threads get the shared variable initializations
  __syncthreads();

  for (int z=blockIdx.x; z < origSize.z; z += gridDim.x) {
    for (int y1=0; y1 < origSize.y; y1 += blockDim.y) {
      int y = y1 + threadIdx.y;
      for (int x1=0; x1 < origSize.x; x1 += blockDim.x) {
        int x = x1 + threadIdx.x;

        if (x < origSize.x && y < origSize.y) {

          float valuef = data[x + dataSize.x*(y + z*dataSize.y)];
          int value = ByteInputData::floatToByte(valuef);
          int origValue = orig[x + origSize.x*(y + z*origSize.y)];

          // printf("%d,%d,%d %d %d\n", z, y, x, origValue, value);

          int diff = abs(value - origValue);
          sumDiff += diff;
          sumDiffSquared += diff*diff;
        }

      }
    }
  }

  // add this thread's results to the total for the thread block
  atomicAdd(&sumDiffShared, sumDiff);
  atomicAdd(&sumDiffSquaredShared, sumDiffSquared);

  // sync so the primary thread gets all the atomicAdds
  __syncthreads();

  // add this thread block's results to errSums[]
  if (threadIdx.x == 0 && threadIdx.y == 0) {
    // FYI, atomicAdd(float) requires CUDA compute level 2.x or higher
    atomicAdd(&errSums[0], (float)sumDiffShared);
    atomicAdd(&errSums[1], (float)sumDiffSquaredShared);
  }

  // overhead from doing atomic updates (with 4 thread blocks of 32x32 threads)
  // without atomic update: 6.501ms
  // with atomic update: 6.742ms
}




void computeErrorRatesAfterDequantGPU
(float *data_dev, scu_wavelet::int3 size, float *tempData_dev,
 const WaveletCompressionParam &param,
 const unsigned char *inputData_dev, ErrorAccumulator &errAccum) {
  
  if (inputData_dev == NULL) {
    fprintf(stderr, "Current implementation can only compute error rates if input data is unsigned bytes (0..255)\n");
    return;
  }

  // inverse transform
  CudaTimer waveletTimer("Inverse wavelet"),
    transposeTimer("Inverse transpose"),
    computeErrTimer("Compute error");

  // printDeviceArray(data_dev, size, "before inverse transform");
  waveletTransformGPU(data_dev, tempData_dev, size, param, true,
                      &waveletTimer, &transposeTimer);
  // printDeviceArray(data_dev, size, "after inverse transform");

  float *errSums_dev;
  computeErrTimer.start();
  CUCHECK(cudaMalloc((void**)&errSums_dev, sizeof(float) * 2));
  CUCHECK(cudaMemset(errSums_dev, 0, sizeof(float) * 2));

  // on W530 laptop, around.cube
  // with 1 thread or (param.originalSize.z) threads, this takes 20ms
  // 2 thread blocks: 11ms
  // 4 thread blocks: 6ms
  // 8 thread blocks: 26ms
  computeErrorKernel<<<4,dim3(32,32)>>>
    (errSums_dev, data_dev, inputData_dev, size, param.originalSize);

  float errSums[2];
  CUCHECK(cudaMemcpy(&errSums, errSums_dev, sizeof(float) * 2,
                     cudaMemcpyDeviceToHost));
  CUCHECK(cudaFree(errSums_dev));
  computeErrTimer.end();

  CUCHECK(cudaThreadSynchronize());
  waveletTimer.print();
  transposeTimer.print();
  computeErrTimer.print();
  // printf("sumDiff = %f, sumDiffSquared = %f\n", errSums[0], errSums[1]);
  
  // largest unsigned char value
  errAccum.setMaxPossible(255);
  errAccum.setSumDiff(errSums[0]);
  errAccum.setSumDiffSquared(errSums[1]);
  errAccum.setCount(param.originalSize.count());

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
  timer.sync();
  timer.print();

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
