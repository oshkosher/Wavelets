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

bool compressFile(const char *inputFile, const char *outputFile, Options &opt);
bool decompressFile(const char *inputFile, const char *outputFile,
		    Options &opt);

bool check1WaveletTransform(const float *data_dev, const CubeFloat *dataBefore,
                            const WaveletCompressionParam &param);

bool check2Sorted(const float *data_dev, const float *orig_data_dev,
                  int dataCount);

bool check3Quantized(const int *quantizedData_dev, int dataCount,
                     const float *unquantizedData_dev,
                     const WaveletCompressionParam &param,
                     float thresholdValue, float max);


void __global__ quantUniformKernel(float *data, int count,
                                   int bits, float threshold, float max);
void __global__ quantLogKernel(float *data, int count,
                               int bits, float threshold, float max);

inline bool isClose(float a, float b) {
  float diff = fabsf(a-b);
  return ((a==0) ? diff : diff/a) < 0.00001;
}

inline bool isClose(int a, int b) {
  return abs(a-b) <= 1;
}

// return the number of mismatches
int compareArrays(const float *a, const float *b, int count);
int compareArrays(const int *a, const int *b, int count);

// debug output
void printArray(float *array, int width, int height, int depth,
                const char *name);

// struct AbsFunctor : public thrust::unary_function<float,float> {
struct AbsFunctor {
  __host__ __device__ float operator() (float x) const {
    return fabsf(x);
  }
};

struct QuantUniformFunctor {
  QuantUniform uni;

  __host__ __device__ QuantUniformFunctor(int bits, float threshold,
                                          float maxVal) {
    uni.init(bits, threshold, maxVal);
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
  double firstStartTime, startTime, elapsed;
  CudaTimer copyToGPUTimer("Copy to GPU"), waveletTimer("Wavelet transform"),
    transposeTimer("Transpose"),
    sortTimer("Sort"), quantizeTimer("Quantize"), 
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
  float *data1_dev, *data2_dev;
  size_t dataCount = data.count();
  size_t dataCountBytes = dataCount * sizeof(float);
  startTime = NixTimer::time();
  CUCHECK(cudaMalloc((void**)&data1_dev, dataCountBytes));
  CUCHECK(cudaMalloc((void**)&data2_dev, dataCountBytes));
  if (!QUIET) {
    printf("CUDA initialize: %.2f ms\n", 1000*(NixTimer::time() - startTime));
    fflush(stdout);
  }

  // copy to GPU
  copyToGPUTimer.start();
  CUCHECK(cudaMemcpy(data1_dev, data.data(), dataCountBytes,
                     cudaMemcpyHostToDevice));
  copyToGPUTimer.end();

  // Wavelet transform
  scu_wavelet::int3 deviceSize = data.size;
  haar_3d_cuda(data1_dev, data2_dev, deviceSize, param.transformSteps,
               false, param.isWaveletTransposeStandard, &waveletTimer,
               &transposeTimer);

  // this is needed to get the final timer value
  CUCHECK(cudaThreadSynchronize());

  // just for testing--copy the data back to the CPU to see if it's correct
  check1WaveletTransform(data1_dev, &data, param);

  // sort absolute values

  // copy the data to data2_dev to be sorted
  CUCHECK(cudaMemcpy(data2_dev, data1_dev, dataCountBytes,
                     cudaMemcpyDeviceToDevice));

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

  AbsFunctor absFunctor;
  sortTimer.start();
  thrust::transform(data1_start, data1_end, data2_start, absFunctor);
  thrust::sort(data2_start, data2_start + dataCount);
  sortTimer.end();
  CUCHECK(cudaThreadSynchronize());

  // check2Sorted(data2_dev, data1_dev, dataCount);

  // data1_dev now contains the data, and data2_dev contains the sorted
  // absolute values of the data

  // fetch percentage threshold values from GPU

  // send back the min, minimum nonzero, and maximum

  float minNonzero, max;
  CUCHECK(cudaMemcpy(&max, data2_dev+dataCount-1, sizeof(float),
                     cudaMemcpyDeviceToHost));


  // if you omit thrust::less<float>(), it defaults to an integer comparison,
  // and looks for the first value greater than or equal to 1
  size_t minNonzeroOffset = 
    thrust::upper_bound(data2_start, data2_end, 0, thrust::less<float>())
    - data2_start;
  CUCHECK(cudaMemcpy(&minNonzero, data2_dev + minNonzeroOffset,
                     sizeof(float), cudaMemcpyDeviceToHost));

  // get the threshold value
  if (param.thresholdFraction <= 0) {
    param.thresholdValue = minNonzero;
  } else if (param.thresholdFraction >= 1) {
    param.thresholdValue = max;
  } else {
    int thresholdOffset = (int) (param.thresholdFraction * dataCount);
    // printf("threshold offest = %d\n", thresholdOffset);
    if (thresholdOffset <= minNonzeroOffset) {
      param.thresholdValue = minNonzero;
    } else {
      CUCHECK(cudaMemcpy(&param.thresholdValue, data2_dev + thresholdOffset,
                         sizeof(float), cudaMemcpyDeviceToHost));
    }
  }

  printf("max = %f, minNonzero = %.7g, threshold = %.10g\n", max,
         minNonzero, param.thresholdValue);
  param.maxValue = max;

  // choose threshold and bin count on CPU

  // quantize on GPU
  float *unquantizedCopy = NULL;

  CUCHECK(cudaMalloc((void**)&unquantizedCopy, dataCountBytes));
  CUCHECK(cudaMemcpy(unquantizedCopy, data1_dev, dataCountBytes,
                     cudaMemcpyDeviceToDevice));

  // overwrite the data in data1_dev with quantized integers

  quantizeTimer.start();
  switch (param.quantAlg) {
  case QUANT_ALG_UNIFORM:
    quantUniformKernel<<<1,1024>>>(data1_dev, dataCount, param.binCount,
                                   param.thresholdValue, max);
    break;
  case QUANT_ALG_LOG:
    quantLogKernel<<<1,1024>>>(data1_dev, dataCount, param.binCount,
                               param.thresholdValue, max);
    break;
  case QUANT_ALG_LLOYD:
    fprintf(stderr, "Lloyd's algorithm not integrated yet.\n");
    return false;
  default:
    fprintf(stderr, "Quantization algorithm %d not found.\n",
            (int)param.quantAlg);
    return false;
  }
  quantizeTimer.end();

  if (unquantizedCopy)
    check3Quantized((int*)data1_dev, dataCount, unquantizedCopy,
                    param, param.thresholdValue, max);

  // copy the data back to the CPU
  copyFromGPUTimer.start();
  quantizedData.size = deviceSize;
  quantizedData.allocate();
  CUCHECK(cudaMemcpy(quantizedData.data(), data1_dev, dataCountBytes,
                     cudaMemcpyDeviceToHost));
  copyFromGPUTimer.end();

  // this is needed to get the final timer value
  CUCHECK(cudaThreadSynchronize());

  copyToGPUTimer.print();
  waveletTimer.print();
  transposeTimer.print();
  sortTimer.print();
  quantizeTimer.print();
  copyFromGPUTimer.print();

  // compute error rates
  if (opt.doComputeError) {
    float meanSqErr, peakSNR;
    // quantizedData.print("quantized before computing error rates");
    computeErrorRates(&quantizedData, param, &inputData, &meanSqErr, &peakSNR);
    printf("Mean squared error: %.3f, peak SNR: %.3f\n", meanSqErr, peakSNR);
  }

  // Huffman code and output

  // compute the quantization bins on the GPU
  // uniform and log: just need the max
  // lloyd: initialize with max, then compute on GPU

  // write the quantized data to a file

  elapsed = NixTimer::time() - firstStartTime;
  printf("Total: %.2f ms\n", elapsed*1000);

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


void __global__ quantUniformKernel(float *data, int count,
                                   int binCount, float threshold, float max) {
  __shared__ QuantUniform quanter;
  int *dataInt = (int*) data;

  if (threadIdx.x == 0) {
    quanter.init(binCount, threshold, max);
    /*
    printf("GPU max = %f, %f\n", max, (double)max);
    printf("GPU bits=%d, thresh=%f, max=%f\n", bits, threshold, max);
    printf("GPU %f->%d\n", 0.0, quanter.quant(0));
    printf("GPU %f->%d\n", 1.0, quanter.quant(1));
    printf("GPU %f->%d\n", .1, quanter.quant(.1));
    printf("GPU %f->%d\n", -.5, quanter.quant(-.5));
    printf("GPU %f->%d\n", 3.0, quanter.quant(3));
    */
  }

  __syncthreads();
  
  int idx = threadIdx.x;
  while (idx < count) {
    dataInt[idx] = quanter.quant(data[idx]);
    idx += blockDim.x;
  }
}


void __global__ quantLogKernel(float *data, int count,
                               int binCount, float threshold, float max) {
  __shared__ QuantLog quanter;
  int *dataInt = (int*) data;

  if (threadIdx.x == 0) {
    // printf("quantLogKernel %d %.10f %.10f\n", binCount, threshold, max);
    quanter.init(binCount, threshold, max);
    // printf("GPU %.10g->%d\n", -0.009803920984, quanter.quant(-0.009803920984));

    /*
    printf("GPU max = %f\n", max);
    printf("GPU bits=%d, thresh=%f, max=%f\n", binCount, threshold, max);
    printf("GPU %f->%d\n", 0.0, quanter.quant(0));
    printf("GPU %f->%d\n", 1.0, quanter.quant(1));
    printf("GPU %f->%d\n", .1, quanter.quant(.1));
    printf("GPU %f->%d\n", -.5, quanter.quant(-.5));
    printf("GPU %f->%d\n", 3.0, quanter.quant(3));
    */
  }

  __syncthreads();
  
  int idx = threadIdx.x;
  while (idx < count) {
    dataInt[idx] = quanter.quant(data[idx]);
    idx += blockDim.x;
  }
}


// return the number of mismatches
int compareArrays(const float *a, const float *b, int count) {
  int mismatches = 0;
  for (int i=0; i < count; i++) {
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


bool check1WaveletTransform(const float *data_dev, const CubeFloat *dataBefore,
                            const WaveletCompressionParam &param) {
  
  // transform the data on the CPU
  CubeFloat result;
  result.size = dataBefore->size;
  result.allocate();
  result.copyFrom(*dataBefore);
  int count = dataBefore->count();

  /*
  printf("Before.\n");
  for (int i=0; i < dataBefore->size.x; i++) {
    printf("%d)\t%.4g\n", i, result.get(i,0,0));
  }
  */

  // printArray(result.data(), result.size.x, result.size.y, result.size.z, "Before");
             
  
  // writeDataFile("check1a.data", data, width, height, true);
  haar_3d(&result, param.transformSteps, false,
          param.isWaveletTransposeStandard);
  // writeDataFile("check1b.data", data, width, height, true);

  // copy the data from the GPU to the CPU
  float *dataTransformedByGPU = new float[count];
  CUCHECK(cudaMemcpy(dataTransformedByGPU, data_dev, sizeof(float) * count,
                     cudaMemcpyDeviceToHost));
  // writeDataFile("check1c.data", dataTransformedByGPU, width, height, true);

  // printArray(result.data(), result.size.x, result.size.y, result.size.z, "After CPU");
  // printArray(dataTransformedByGPU, result.size.x, result.size.y, result.size.z, "After GPU");

  /*
  printf("CPU.\n");
  for (int i=0; i < dataBefore->size.x; i++) {
    printf("%d)\t%.4g\n", i, result.get(i,0,0));
  }

  printf("GPU.\n");
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
                     float thresholdValue, float max) {

  float *unquantized = new float[count];
  CUCHECK(cudaMemcpy(unquantized, unquantizedData_dev, count * sizeof(float),
                     cudaMemcpyDeviceToHost));

  int *gpuQuant = new int[count];
  CUCHECK(cudaMemcpy(gpuQuant, quantizedData_dev, count * sizeof(int),
                     cudaMemcpyDeviceToHost));
  
  Quantizer *quanter = createQuantizer(param);
  int *cpuQuant = new int[count];
  quanter->quantizeRow(unquantized, cpuQuant, count);

  int mismatches = compareArrays(cpuQuant, gpuQuant, count);
  int notEqual = 0;
  for (int i=0; i < count; i++) {
    if (cpuQuant[i] != gpuQuant[i]) {
      notEqual++;
      // printf("[%d] %.10g  %d %d\n", i, unquantized[i], cpuQuant[i], gpuQuant[i]);
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
