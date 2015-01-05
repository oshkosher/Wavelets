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

bool check1WaveletTransform(const float *data_dev, const float *dataBefore,
                            int width, int height, const Options &opt);

bool check2Sorted(const float *data_dev, const float *orig_data_dev,
                  int dataCount);

bool check3Quantized(const int *quantizedData_dev, int dataCount,
                     const float *unquantizedData, const Options &opt,
                     float thresholdValue, float max);


void __global__ quantUniformKernel(float *data, int count,
                                   int bits, float threshold, float max);
void __global__ quantLogKernel(float *data, int count,
                               int bits, float threshold, float max);

inline bool isClose(float a, float b) {
  float diff = fabsf(a-b);
  return ((a==0) ? diff : diff/a) < 0.00001;
}

// return the number of mismatches
int compareArrays(const float *a, const float *b, int count);

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

  if (argc - nextArg != 2) printHelp();

  const char *inputFile = argv[nextArg++];
  const char *outputFile = argv[nextArg++];
  bool result;

  if (opt.doCompress) {
    result = compressFile(inputFile, outputFile, opt);
  } else {
    result = decompressFile(inputFile, outputFile, opt);
  }
  if (result == false) return 1;

  return 0;
}

bool compressFile(const char *inputFile, const char *outputFile,
                  Options &opt) {

  // read the data file
  float *data;
  int width, height;
  double firstStartTime, startTime, elapsed;
  CudaTimer copyToGPUTimer("Copy to GPU"), waveletTimer("wavelet transform"),
    sortTimer("Sort"), quantizeTimer("Quantize"), 
    copyFromGPUTimer("Copy from GPU");

  firstStartTime = startTime = NixTimer::time();

  if (!readDataFile(inputFile, &data, &width, &height)) return 1;
  elapsed = NixTimer::time() - startTime;
  printf("Read %dx%d data file: %.2f ms\n", width, height, elapsed * 1000);

  // pad the data to make it a square power of two, if necessary
  // XXX implement this on the GPU
  startTime = NixTimer::time();
  int longSide = width > height ? width : height;
  int size = dwt_padded_length(longSide, 0, true);
  if (width != size || height != size) {
    startTime = NixTimer::time();
    float *paddedData = dwt_pad_2d(height, width, width, data,
                                   size, size, size,
                                   NULL, REFLECT);
    elapsed = NixTimer::time() - startTime;
    printf("Pad data: %.2f ms\n", elapsed*1000);
    delete[] data;
    data = paddedData;
    width = height = size;
  }
  elapsed = NixTimer::time() - startTime;
  printf("Pad data to %dx%d: %.2f ms\n", width, height, elapsed * 1000);

  // the data is now a square, size width*height

  // adjust the number of wavelet steps in case the user requested too many
  int maxWaveletSteps = dwtMaximumSteps(size);
  if (opt.waveletSteps > maxWaveletSteps)
    opt.waveletSteps = maxWaveletSteps;

  // send the data to the GPU
  // XXX use a stream for asynchronous operation
  float *data1_dev, *data2_dev;
  size_t dataCount = (size_t) size * size;
  size_t dataCountBytes = dataCount * sizeof(float);
  CUCHECK(cudaMalloc((void**)&data1_dev, dataCountBytes));
  CUCHECK(cudaMalloc((void**)&data2_dev, dataCountBytes));

  // copy to GPU
  copyToGPUTimer.start();
  CUCHECK(cudaMemcpy(data1_dev, data, dataCountBytes, cudaMemcpyHostToDevice));
  copyToGPUTimer.end();

  // Wavelet transform
  waveletTimer.start();
  wavelet_transform_gpu(data1_dev, data2_dev, size, opt.waveletSteps,
                        false, opt.isWaveletTransposeStandard);
  waveletTimer.end();
  
  // data in data1_dev is modified in-place

  // just for testing--copy the data back to the CPU to see if it's correct
  // check1WaveletTransform(data1_dev, data, size, size, opt);

  // copy it to data2_dev to be sorted

  // sort the absolute values on the GPU
  thrust::device_ptr<float>
    data1_start(data1_dev),
    data1_end(data1_dev + dataCount),
    data2_start(data2_dev),
    data2_end(data2_dev + dataCount);

  // CUCHECK(cudaMemcpy(data2_dev, data1_dev, dataCountBytes,
  // cudaMemcpyDeviceToDevice));

  /* thrust::sort works with pointers on the device (like data2_dev), 
     as long as you give it "thrust::device" as the first argument.
     thrust::transform does not. It requires thrust::device_ptr objects
     (like data1_start). thrust::sort will work with those too.
  */

  // make a copy of the data, applying abs() to each value, into data2,
  // then sort that copy
  sortTimer.start();
  AbsFunctor absFunctor;
  thrust::transform(data1_start, data1_end, data2_start,
                    // data1_dev, data1_dev + dataCount, data2_dev,
                    absFunctor);
  thrust::sort(data2_start, data2_start + dataCount);
  sortTimer.end();

  // XXX try out sorting with a transform_iterator

  // data1_dev now contains the data, and data2_dev contains the sorted
  // absolute values of the data

  // check2Sorted(data2_dev, data1_dev, dataCount);

  // send back the min, minimum nonzero, and maximum
  float min, minNonzero, max, thresholdValue;
  CUCHECK(cudaMemcpy(&min, data2_dev, sizeof(float), cudaMemcpyDeviceToHost));
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
  if (opt.thresholdFraction <= 0) {
    thresholdValue = min;
  } else if (opt.thresholdFraction >= 1) {
    thresholdValue = max;
  } else {
    int thresholdOffset = (int) (opt.thresholdFraction * dataCount);
    // printf("threshold offest = %d\n", thresholdOffset);
    if (thresholdOffset <= minNonzeroOffset) {
      thresholdValue = minNonzero;
    } else {
      CUCHECK(cudaMemcpy(&thresholdValue, data2_dev + thresholdOffset,
                         sizeof(float), cudaMemcpyDeviceToHost));
    }
  }

  printf("min = %f, max = %f, minNonzero = %.7g, threshold = %.7g\n", min, max,
         minNonzero, thresholdValue);

  std::vector<float> quantBinBoundaries;
  std::vector<float> quantBinValues;

  // compute the quantization bins on the GPU
  // uniform and log: just need the max
  // lloyd: initialize with max, then compute on GPU
  // count: probably no faster to multithread

  float *unquantizedCopy = NULL;
  /*
  CUCHECK(cudaMalloc((void**)&unquantizedCopy, dataCountBytes));
  CUCHECK(cudaMemcpy(unquantizedCopy, data1_dev, dataCountBytes,
                     cudaMemcpyDeviceToDevice));
  */

  quantizeTimer.start();
  switch (opt.quantizeAlgorithm) {
  case QUANT_ALG_UNIFORM:
    quantUniformKernel<<<1,1024>>>(data1_dev, dataCount,
                                   opt.quantizeBits, thresholdValue, max);
    break;
  case QUANT_ALG_LOG:
    quantLogKernel<<<1,1024>>>(data1_dev, dataCount,
                               opt.quantizeBits, thresholdValue, max);
    break;
  case QUANT_ALG_COUNT:
    fprintf(stderr, "Count algorithm not integrated yet.\n");
    /*
    quant_count_init_sorted_cpu(size*size, sortedAbsData, opt.quantizeBits,
				threshold, quantBinBoundaries, quantBinValues);
    quant_boundaries_array(quantBinBoundaries, size*size, data);
    */
    break;
  case QUANT_ALG_LLOYD:
    fprintf(stderr, "Lloyd's algorithm not integrated yet.\n");
    return false;
  default:
    fprintf(stderr, "Quantization algorithm %d not found.\n",
            (int)opt.quantizeAlgorithm);
    return false;
  }
  quantizeTimer.end();

  if (unquantizedCopy)
    check3Quantized((int*)data1_dev, dataCount, unquantizedCopy,
                    opt, thresholdValue, max);

  // XXX perform further compressions of the data?

  // copy the data back to the CPU
  copyFromGPUTimer.start();
  CUCHECK(cudaMemcpy(data, data1_dev, dataCountBytes, cudaMemcpyDeviceToHost));
  copyFromGPUTimer.end();

  // this is needed to get the final timer value
  CUCHECK(cudaThreadSynchronize());

  copyToGPUTimer.print();
  waveletTimer.print();
  sortTimer.print();
  quantizeTimer.print();
  copyFromGPUTimer.print();

  // write the quantized data to a file
  FileData fileData(opt, NULL, (int*)data, size, size);
  fileData.threshold = thresholdValue;
  if (opt.quantizeAlgorithm == QUANT_ALG_UNIFORM ||
      opt.quantizeAlgorithm == QUANT_ALG_LOG) {
    fileData.quantMaxVal = max;
  } else {
    fileData.quantBinBoundaries = quantBinBoundaries;
    fileData.quantBinValues = quantBinValues;
  }

  startTime = NixTimer::time();
  if (!writeQuantData(outputFile, fileData)) return false;
  elapsed = NixTimer::time() - startTime;
  printf("Write data file: %.2f ms\n", elapsed*1000);

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


bool wavelet_transform_gpu(float *data_dev, float *tmp_data_dev,
                           int size, int stepCount, bool isInverse,
                           bool isStandardTranspose, cudaStream_t stream) {

  // XXX isStandardTranspose not implemented. Always does nonstandard.
  if (isStandardTranspose)
    printf("FYI, isStandardTranspose not implemented yet.\n");


  int tileSize = bestHaarGPUTileSize();
  size_t sharedMemSize = tileSize * (tileSize+1)
    * 2 * sizeof(float);

  if (!isInverse) {

    int transformLength = size;

    for (int i=0; i < stepCount; i++) {

      dim3 gridDim((transformLength - 1) / (tileSize*2) + 1,
                   (transformLength - 1) / (tileSize) + 1);
      dim3 blockDim(tileSize, tileSize);
      // printf("grid: %dx%d, block: %dx%d\n", (int)gridDim.x, (int)gridDim.y,
      // (int)blockDim.x, (int)blockDim.y);

      // do the wavelet transform on rows  data->tmp
      haar_transpose_2d_kernel
        <<<gridDim, blockDim, sharedMemSize, stream>>>
        (size, transformLength, data_dev, tmp_data_dev, tileSize);
      // printf("transpose kernel rows %d\n", i);

      // do the wavelet transform on columns  tmp->data
      haar_transpose_2d_kernel
        <<<gridDim, blockDim, sharedMemSize, stream>>>
        (size, transformLength, tmp_data_dev, data_dev, tileSize);
      // printf("transpose kernel columns %d\n", i);

      transformLength >>= 1;
    }

  } else {

    // isInverse==true
    int transformLength = size >> (stepCount - 1);

    for (int i=0; i < stepCount; i++) {

      dim3 gridDim((transformLength - 1) / (tileSize*2) + 1,
                   (transformLength - 1) / (tileSize) + 1);
      dim3 blockDim(tileSize, tileSize);
    
      // transform columns and transpose  data->tmp
      haar_inv_transpose_2d_kernel
        <<<gridDim, blockDim, sharedMemSize, stream>>>
        (size, transformLength, data_dev, tmp_data_dev, tileSize);
    
      // transform rows and transpose  tmp->data
      haar_inv_transpose_2d_kernel
        <<<gridDim, blockDim, sharedMemSize, stream>>>
        (size, transformLength, tmp_data_dev, data_dev, tileSize);
    
      transformLength <<= 1;
    }
  }
  
  return true;
}


void __global__ quantUniformKernel(float *data, int count,
                                   int bits, float threshold, float max) {
  __shared__ QuantUniform quanter;
  int *dataInt = (int*) data;

  if (threadIdx.x == 0) {
    quanter.init(bits, threshold, max);
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
                               int bits, float threshold, float max) {
  __shared__ QuantLog quanter;
  int *dataInt = (int*) data;

  if (threadIdx.x == 0) {
    quanter.init(bits, threshold, max);
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


// return the number of mismatches
int compareArrays(const float *a, const float *b, int count) {
  int mismatches = 0;
  for (int i=0; i < count; i++) {
    if (!isClose(a[i], b[i])) mismatches++;
  }
  return mismatches;
}


bool check1WaveletTransform(const float *data_dev, const float *dataBefore,
                            int width, int height, const Options &opt) {
  
  // transform the data on the CPU
  float *data = new float[width*height];
  memcpy(data, dataBefore, sizeof(float) * width * height);
  writeDataFile("check1a.data", data, width, height, true);
  haar_2d(width, data, false, opt.waveletSteps);
  writeDataFile("check1b.data", data, width, height, true);

  // copy the data from the GPU to the CPU
  float *dataTransformedByGPU = new float[width*height];
  CUCHECK(cudaMemcpy(dataTransformedByGPU, data_dev,
                     sizeof(float) * width * height, cudaMemcpyDeviceToHost));
  writeDataFile("check1c.data", dataTransformedByGPU, width, height, true);

  // compare the data
  int mismatches = compareArrays(data, dataTransformedByGPU, width*height);

  printf("Check 1: wavelet transform %d / %d mismatches\n",
         mismatches, width*height);

  delete[] data;
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
                     const float *unquantizedData_dev, const Options &opt,
                     float thresholdValue, float max) {

  int *quantizedData = new int[count];
  CUCHECK(cudaMemcpy(quantizedData, quantizedData_dev, count * sizeof(int),
                     cudaMemcpyDeviceToHost));
  float *unquantizedData = new float[count];
  CUCHECK(cudaMemcpy(unquantizedData, unquantizedData_dev,
                     count * sizeof(float),
                     cudaMemcpyDeviceToHost));
  /*
  for (int i=0; i < count; i++) {
    printf("%d\t%f\n", quantizedData[i], unquantizedData[i]);
  }
  */
  delete[] quantizedData;
  delete[] unquantizedData;

  return true;
}

  
  


