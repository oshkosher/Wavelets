#include "cuda.h"
#include <thrust/execution_policy.h>
#include <thrust/functional.h>
#include <thrust/transform.h>
#include <thrust/sort.h>
#include "test_compress_common.h"
#include "dquant_log_cpu.h"
#include "dquant_unif_cpu.h"
#include "dwt_cpu.h"
#include "nixtimer.h"
#include "quant_count.h"
#include "quant_log_cpu.h"
#include "quant_unif_cpu.h"
#include "thresh_cpu.h"
#include "cucheck.h"

bool compressFile(const char *inputFile, const char *outputFile, Options &opt);
bool decompressFile(const char *inputFile, const char *outputFile,
		    Options &opt);


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

struct AbsFunctor : public thrust::unary_function<float,float> {
  __device__ float operator() (float x) const {
    return fabsf(x);
  }
};

bool compressFile(const char *inputFile, const char *outputFile,
                  Options &opt) {

  // read the data file
  float *data;
  int width, height;
  double firstStartTime, startTime, elapsed;

  firstStartTime = startTime = NixTimer::time();

  if (!readDataFile(inputFile, &data, &width, &height)) return 1;
  elapsed = NixTimer::time() - startTime;
  printf("Read %dx%d data file: %.2f ms\n", width, height, elapsed * 1000);

  // pad the data to make it a square power of two, if necessary
  // XXX implement this on the GPU
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
  CUCHECK(cudaMemcpy(data1_dev, data, dataCountBytes, cudaMemcpyHostToDevice));

  float waveletMs = 0;
  // wavelet_transform_gpu(data1_dev, data2_dev, size, size, opt.waveletSteps, opt.isWaveletTransposeStandard, waveletMs);
                        

  // data in data1_dev is modified in-place

  // perform the wavelet transform
  // float waveletMs = haar_2d(size, data, false, opt.waveletSteps);
  printf("Wavelet transform: %.2f ms\n", waveletMs);

  // sort the absolute values on the GPU
  /*thrust::device_ptr<float>
    data1_start(data1_dev),
    data1_end(data1_dev + dataCount),
    data2_start(data2_dev),
    data2_end(data2_dev + dataCount);*/
  AbsFunctor absFunctor;
  thrust::transform(thrust::device, data1_dev, data1_dev + dataCount,
                    data2_dev, absFunctor);
  thrust::sort(thrust::device, data2_dev, data2_dev + dataCount);

  // send back the max

  // find the threshold value on the GPU (single-thread kernel?)
  /*
  float maxVal, *sortedAbsData = NULL;
  startTime = NixTimer::time();
  float threshold = thresh_cpu(size, data, opt.thresholdFraction,
                               &maxVal, &sortedAbsData);
  elapsed = NixTimer::time() - startTime;
  */
  // printf("threshold = %g: %.2f ms\n", threshold, elapsed*1000);

  std::vector<float> quantBinBoundaries;
  std::vector<float> quantBinValues;

  // compute the quantization bins on the GPU
  // uniform and log: just need the max
  // lloyd: initialize with max, then compute on GPU
  // count: probably no faster to multithread

  // threshold/quantize the data

  // quantize the data

  /*
  startTime = NixTimer::time();
  switch (opt.quantizeAlgorithm) {
  case QUANT_ALG_UNIFORM:
    quant_unif_cpu(size, data, opt.quantizeBits, threshold, maxVal);
    break;
  case QUANT_ALG_LOG:
    quant_log_cpu(size, data, opt.quantizeBits, threshold, maxVal);
    break;
  case QUANT_ALG_COUNT:
    quant_count_init_sorted_cpu(size*size, sortedAbsData, opt.quantizeBits,
				threshold, quantBinBoundaries, quantBinValues);
    quant_boundaries_array(quantBinBoundaries, size*size, data);
    break;
  case QUANT_ALG_LLOYD:
    fprintf(stderr, "Lloyd's algorithm not integrated yet.\n");
    return false;
  default:
    fprintf(stderr, "Quantization algorithm %d not found.\n",
            (int)opt.quantizeAlgorithm);
    return false;
  }
  elapsed = NixTimer::time() - startTime;
  printf("Apply threshold / quantize: %.2f ms\n", elapsed*1000);

  delete[] sortedAbsData;

  // perform further compressions of the data?

  // copy the data back to the CPU

  // write the quantized data to a file
  FileData fileData(opt, data, size, size);
  fileData.threshold = threshold;
  if (opt.quantizeAlgorithm == QUANT_ALG_UNIFORM ||
      opt.quantizeAlgorithm == QUANT_ALG_LOG) {
    fileData.quantMaxVal = maxVal;
  } else {
    fileData.quantBinBoundaries = quantBinBoundaries;
    fileData.quantBinValues = quantBinValues;
  }

  startTime = NixTimer::time();

  // if (!writeQuantDataSimple(outputFile, fileData)) return false;
  // if (!writeQuantDataParamStrings(outputFile, fileData)) return false;
  if (!writeQuantDataProtoBuf(outputFile, fileData)) return false;
    
  elapsed = NixTimer::time() - startTime;
  printf("Write data file: %.2f ms\n", elapsed*1000);

  elapsed = NixTimer::time() - firstStartTime;
  printf("Total: %.2f ms\n", elapsed*1000);
  */
  

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

  // if (!readQuantDataSimple(inputFile, f)) return false;
  // if (!readQuantDataParamStrings(inputFile, f)) return false;
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
