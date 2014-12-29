#include "test_compress_common.h"
#include "dwt_cpu.h"
#include "nixtimer.h"
#include "thresh_cpu.h"
#include "quant.h"
#include "lloyds.h"

bool compressFile(const char *inputFile, const char *outputFile, Options &opt);
bool decompressFile(const char *inputFile, const char *outputFile,
		    Options &opt);
void computeLloydQuantization(const float *inputData, int count, 
			      float minVal, float maxVal, int bits,
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
  int count = size * size;
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
  if (opt.waveletSteps > maxWaveletSteps || opt.waveletSteps < 0)
    opt.waveletSteps = maxWaveletSteps;

  // perform the wavelet transform
  if (opt.waveletSteps > 0) {
    float waveletMs = haar_2d(size, data, false, opt.waveletSteps);
    printf("Wavelet transform (%d steps): %.2f ms\n", 
           opt.waveletSteps, waveletMs);
  }

  // save the intermediate data to a file before quantizing
  if (opt.saveBeforeQuantizingFilename != "") {
    const char *filename = opt.saveBeforeQuantizingFilename.c_str();
    if (!writeDataFile(filename, data, width, height)) {
                       
      printf("Failed to write intermediate data file \"%s\".\n", filename);
    } else {
      printf("Write intermediate data file \"%s\"\n", filename);
    }
  }

  // find the threshold value by sorting
  // XXX use quickselect to speed up
  float maxVal, minVal, maxAbsVal, *sortedAbsData;
  int nonzeroCount;
  // float *sortedAbsData = NULL;
  startTime = NixTimer::time();
  float threshold = thresh_cpu(count, data, opt.thresholdFraction,
			       &nonzeroCount, &maxVal, &minVal,
			       &sortedAbsData);
			       
  maxAbsVal = sortedAbsData[count-1];
  float *nonzeroData = sortedAbsData + count - nonzeroCount;

  elapsed = NixTimer::time() - startTime;
  printf("Compute threshold = %g, min = %g, max = %g: %.2f ms\n",
	 threshold, minVal, maxVal, elapsed*1000);

  std::vector<float> quantBinBoundaries;
  std::vector<float> quantBinValues;
  int *quantizedData = NULL;

  // quantize the data
  startTime = NixTimer::time();
  switch (opt.quantizeAlgorithm) {

  case QUANT_ALG_UNIFORM:
    // quant_unif_cpu(size, data, opt.quantizeBits, threshold, maxVal);
    {
      QuantUniform qunif(opt.quantizeBits, threshold, maxAbsVal);
      QuantizationLooper<QuantUniform> qloop(&qunif);
      quantizedData = new int[size*size];
      qloop.quantize(size*size, data, quantizedData, true);
      printf("Quantization error: %g\n", qloop.getError());
    }
    break;

  case QUANT_ALG_LOG:
    // quant_log_cpu(size, data, opt.quantizeBits, threshold, maxVal);
    {
      QuantLog qlog(opt.quantizeBits, threshold, maxAbsVal);
      QuantizationLooper<QuantLog> qloop(&qlog);
      quantizedData = new int[size*size];
      qloop.quantize(size*size, data, quantizedData, true);
      printf("Quantization error: %g\n", qloop.getError());
    }
    break;

  case QUANT_ALG_COUNT:
    /*
    quant_count_init_sorted_cpu(size*size, sortedAbsData, opt.quantizeBits,
				threshold, quantBinBoundaries, quantBinValues);
    quant_boundaries_array(quantBinBoundaries, size*size, data);
    */
    {
      QuantCodebook qcb;
      qcb.initCountBins(count, data, opt.quantizeBits, threshold);
      quantBinBoundaries = qcb.boundaries;
      quantBinValues = qcb.codebook;
      QuantizationLooper<QuantCodebook> qloop(&qcb);
      quantizedData = new int[count];
      qloop.quantize(count, data, quantizedData, true);
      printf("Quantization error: %g\n", qloop.getError());
    }
    break;

  case QUANT_ALG_LLOYD: 
    {
      computeLloydQuantization(nonzeroData, nonzeroCount,
			       threshold, maxAbsVal, opt.quantizeBits,
			       quantBinBoundaries, quantBinValues);
      QuantCodebook qcb(quantBinBoundaries, quantBinValues);
      QuantizationLooper<QuantCodebook> qloop(&qcb);
      quantizedData = new int[count];
      qloop.quantize(count, data, quantizedData, true);
      printf("Quantization error: %g\n", qloop.getError());
      
      delete[] data;
      data = NULL;
    }
    break;

  default:
    fprintf(stderr, "Quantization algorithm %d not found.\n",
            (int)opt.quantizeAlgorithm);
    return false;
  }

  elapsed = NixTimer::time() - startTime;
  printf("Quantize: %.2f ms\n", elapsed*1000);

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
  int *inputData;
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

  inputData = f.intData;
  width = f.width;
  height = f.height;

  assert(inputData != NULL && width > 0 && height > 0);
  assert(width == height);
  int size = width;
  int count = width * height;
  data = new float[count];

  // de-quantize the data
  startTime = NixTimer::time();
  switch (f.quantizeAlgorithm) {
  case QUANT_ALG_UNIFORM:
    // dquant_unif_cpu(size, data, f.quantizeBits, f.threshold, f.quantMaxVal);
    {
      QuantUniform qunif(f.quantizeBits, f.threshold, f.quantMaxVal);
      QuantizationLooper<QuantUniform> qloop(&qunif);
      qloop.dequantize(count, inputData, data);
    }
    break;

  case QUANT_ALG_LOG:
    // dquant_log_cpu(size, data, f.quantizeBits, f.threshold, f.quantMaxVal);
    {
      QuantLog qunif(f.quantizeBits, f.threshold, f.quantMaxVal);
      QuantizationLooper<QuantLog> qloop(&qunif);
      qloop.dequantize(count, inputData, data);
    }
    break;

  case QUANT_ALG_COUNT:
  case QUANT_ALG_LLOYD:
    // dequant_codebook_array(f.quantBinValues, size*size, data);
    {
      QuantCodebook qcb;
      qcb.init(f.quantBinBoundaries, f.quantBinValues);
      QuantizationLooper<QuantCodebook> qloop(&qcb);
      qloop.dequantize(count, inputData, data);
    }
    
    break;

  default:
    fprintf(stderr, "Quantization algorithm %d not found.\n",
            (int)f.quantizeAlgorithm);
    return false;
  }

  delete[] inputData;
  f.intData = NULL;

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


/*
  Distribute N bins like this:

                                    1 bin
                                      |
    N/2-1 bins    1 bin   N/2-1 bins  v
  ---negatives--+-------+--positives--x
                ^   0   ^             ^
                |       |             |
          -thresh       +thresh      max

  The negative and positive bins will be mirrored.

  For example, if N = 8, thresh=1, and max=10, one possible
  set of positive thresholds is: 3, 6
  And the total set of thresholds will be:
    -6 -3 -1 1 3 6 10
  Possible codebook values:
   -8 -5 -2 0 2 5 8 10
*/
void computeLloydQuantization(const float *inputData, int count, 
			      float minVal, float maxVal, int bits,
			      std::vector<float> &quantBinBoundaries,
			      std::vector<float> &quantBinValues) {

  int binCount = (1 << (bits - 1)) - 1;
  quantBinBoundaries.clear();
  quantBinValues.clear();
  float *binBoundaries = new float[binCount-1];
  float *binValues = new float[binCount];

  assert(minVal > 0);
  assert(maxVal > minVal);

  // use log quantization to create an initial codebook
  QuantLog qlog(bits, minVal, maxVal);
  for (int i=0; i < binCount; i++) {
    float dq = qlog.dequant(i+1);
    binValues[i] = dq;
  }

  /*
  // apply the threshold
  float *data = new float[count];
  for (unsigned i=0; i < count; i++) {
    float f = inputData[i];
    if (fabsf(f) <= threshold) f = 0;
    data[i] = f;
  }
  */

  /*
  for (int i=0; i < binCount; i++) {
    printf("InitCB %d: %f\n", i, binValues[i]);
  }
  */


  // fine-tune the codebook and bin boundaries using Lloyd's algorithm.
  // This also applies the quantization to each value, writing the values
  // to quantizedData[]
  float dist, reldist;
  unsigned *quantizedData = new unsigned[count];
  lloyd(inputData, count, binValues, binCount, binBoundaries, dist,
        reldist, quantizedData);

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

  // negative bins
  quantBinValues.push_back(-binValues[binCount-1]);

  for (int i=binCount-2; i >= 0; i--) {
    quantBinBoundaries.push_back(-binBoundaries[i]);
    quantBinValues.push_back(-binValues[i]);
  }

  // zero bin
  quantBinBoundaries.push_back(-minVal);
  quantBinValues.push_back(0);
  quantBinBoundaries.push_back(minVal);

  // positive bins
  for (int i=0; i < binCount-1; i++) {
    quantBinValues.push_back(binValues[i]);
    quantBinBoundaries.push_back(binBoundaries[i]);
  }    
  quantBinValues.push_back(binValues[binCount-1]);

  // top bin
  quantBinBoundaries.push_back(maxVal);
  quantBinValues.push_back(maxVal);

  /*
  for (size_t i = 0; ; i++) {
    printf("Bin %3d. %f\n", (int)i, quantBinValues[i]);
    if (i == quantBinBoundaries.size()) break;
    printf("  %f\n", quantBinBoundaries[i]);
  }
  */

  /*
  for (int i=0; i < count; i++) {
    printf("%g\t%d\n", inputData[i], quantizedData[i]);
  }
  */
  
  delete[] quantizedData;
  delete[] binBoundaries;
  delete[] binValues;
}

  
