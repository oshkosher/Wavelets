#include "test_compress_common.h"
#include "dquant_log_cpu.h"
#include "dquant_unif_cpu.h"
#include "dwt_cpu.h"
#include "nixtimer.h"
#include "quant_count.h"
#include "quant_log_cpu.h"
#include "quant_unif_cpu.h"
#include "thresh_cpu.h"
#include "quant.h"
#include "lloyds.h"

bool compressFile(const char *inputFile, const char *outputFile, Options &opt);
bool decompressFile(const char *inputFile, const char *outputFile,
		    Options &opt);
void applyLloydQuantization(float *inputData, unsigned count, 
                            int *quantizedData, float threshold,
                            float maxAbsVal, int bits,
                            std::vector<float> &quantBinBoundaries,
                            std::vector<float> &quantBinValues);


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

  firstStartTime = startTime = NixTimer::time();

  if (!readDataFile(inputFile, &data, &width, &height)) return 1;
  elapsed = NixTimer::time() - startTime;
  printf("Read %dx%d data file: %.2f ms\n", width, height, elapsed * 1000);

  // pad the data to make it a square power of two, if necessary
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

  // perform the wavelet transform
  float waveletMs = haar_2d(size, data, false, opt.waveletSteps);
  printf("Wavelet transform: %.2f ms\n", waveletMs);

  // find the threshold value
  float maxVal, *sortedAbsData = NULL;
  startTime = NixTimer::time();
  float threshold = thresh_cpu(size, data, opt.thresholdFraction,
                               &maxVal, &sortedAbsData);
  elapsed = NixTimer::time() - startTime;
  printf("threshold = %g: %.2f ms\n", threshold, elapsed*1000);

  std::vector<float> quantBinBoundaries;
  std::vector<float> quantBinValues;
  int *quantizedData = NULL;

  // quantize the data
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
    quantizedData = new int[size*size];
    applyLloydQuantization(data, size*size, quantizedData, threshold, maxVal,
                           opt.quantizeBits, quantBinBoundaries,
                           quantBinValues);
    delete[] data;
    data = NULL;
    break;

  default:
    fprintf(stderr, "Quantization algorithm %d not found.\n",
            (int)opt.quantizeAlgorithm);
    return false;
  }

  elapsed = NixTimer::time() - startTime;
  printf("Apply threshold / quantize: %.2f ms\n", elapsed*1000);

  delete[] sortedAbsData;

  // write the quantized data to a file
  FileData fileData(opt, data, quantizedData, size, size);
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

  delete[] data;

  return true;
}

bool decompressFile(const char *inputFile, const char *outputFile,
		    Options &opt) {

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
  case QUANT_ALG_LLOYD:
    dequant_codebook_array(f.quantBinValues, size*size, data);
    break;
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

  delete[] data;

  return true;
}


void applyLloydQuantization(float *inputData, unsigned count, 
                            int *quantizedData, float threshold,
                            float maxAbsVal, int bits,
                            std::vector<float> &quantBinBoundaries,
                            std::vector<float> &quantBinValues) {

  int binCount = 1 << bits;
  quantBinBoundaries.clear();
  quantBinValues.clear();
  float *binBoundaries = new float[binCount-1];
  float *binValues = new float[binCount];

  assert(maxAbsVal > 0);
  assert(threshold >= 0);
  assert(threshold < maxAbsVal);

  // use log quantization to create an initial codebook
  QuantLog qlog;
  qlog.init(bits, threshold, maxAbsVal);
  int half = binCount/2;
  binValues[half-1] = threshold / -2;
  binValues[half] = threshold / 2;
  for (int i=1; i < half; i++) {
    float dq = qlog.dequant(i);
    binValues[half+i] = dq;
    binValues[half-i-1] = -dq;
  }

  // apply the threshold
  for (unsigned i=0; i < count; i++)
    if (inputData[i] <= threshold) inputData[i] = 0;

  /*
  for (int i=0; i < binCount; i++) {
    printf("InitCB %d: %f\n", i, binValues[i]);
  }
  */

  // fine-tune the codebook and bin boundaries using Lloyd's algorithm.
  // This also applies the quantization to each value, writing the values
  // to quantizedData[]
  float dist, reldist;
  lloyd(inputData, count, binValues, binCount, binBoundaries, dist,
        reldist, (unsigned*) quantizedData);
  // printf("dist = %f, reldist = %f\n", dist, reldist);

  // sanity-check
  for (int i=0; i < binCount-1; i++) {
    if (binValues[i] > binValues[i+1]) {
      printf("ERROR: codebook[%d] > codebook[%d]  (%f > %f)\n",
             i, i+1, binValues[i], binValues[i+1]);
    }

    if (binBoundaries[i] < binValues[i] || binBoundaries[i] > binValues[i+1]) {
      printf("ERROR: partition[%d] (%f) should be between codebook[%d] (%f) and codebook[%d] (%f)\n",
             i, binBoundaries[i], i, binValues[i], i+1, binValues[i+1]);
    }
  }
  
  for (int i=0; i < binCount; i++) {
    quantBinValues.push_back(binValues[i]);
    // printf("Bin %3d. %f\n", i, binValues[i]);
    if (i == binCount-1) break;
    quantBinBoundaries.push_back(binBoundaries[i]);
    // printf("  %f\n", binBoundaries[i]);
  }

  /*
  for (unsigned i=0; i < count; i++) {
    printf("%f\t%d\n", inputData[i], quantizedData[i]);
  }
  */

  delete[] binBoundaries;
  delete[] binValues;
}

  
