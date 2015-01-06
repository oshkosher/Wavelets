#include <google/protobuf/stubs/common.h>
#include "test_compress_common.h"
#include "dwt_cpu.h"
#include "nixtimer.h"
#include "thresh_cpu.h"
#include "quant.h"
#include "lloyds.h"
#include "bit_stream.h"
#include "huffman.h"

bool compressFile(const char *inputFile, const char *outputFile, Options &opt);
bool decompressFile(const char *inputFile, const char *outputFile,
		    Options &opt);
void computeLloydQuantization(const float *inputData, int count, int bits,
			      std::vector<float> &quantBinBoundaries,
			      std::vector<float> &quantBinValues);
void testHuffman(const int data[], int count, int bitCount);

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

  // deallocate static protobuf data
  google::protobuf::ShutdownProtobufLibrary();

  if (result == false) return 1;

  return 0;
}


bool compressFile(const char *inputFile, const char *outputFile,
                  Options &opt) {

  Data2d data, quantizedData;
  double firstStartTime, startTime, elapsed;

  firstStartTime = startTime = NixTimer::time();

  // read the data file
  if (!readDataFile(inputFile, &data.floatData, &data.width, &data.height))
    return 1;
  elapsed = NixTimer::time() - startTime;
  printf("Read %dx%d data file: %.2f ms\n", data.width, data.height,
         elapsed * 1000);

  // pad the data to make it a square power of two, if necessary
  int longSide = data.width > data.height ? data.width : data.height;
  int size = dwt_padded_length(longSide, 0, true);
  int count = size * size;
  if (data.width != size || data.height != size) {
    startTime = NixTimer::time();
    float *paddedData = dwt_pad_2d(data.height, data.width, data.width,
                                   data.floatData,
                                   size, size, size,
                                   NULL, REFLECT);
    elapsed = NixTimer::time() - startTime;
    printf("Pad data: %.2f ms\n", elapsed*1000);
    delete[] data.floatData;
    data.floatData = paddedData;
    data.width = data.height = size;
  }

  const int width = data.width;
  const int height = data.height;

  // adjust the number of wavelet steps in case the user requested too many
  int maxWaveletSteps = dwtMaximumSteps(size);
  if (opt.waveletSteps > maxWaveletSteps || opt.waveletSteps < 0)
    opt.waveletSteps = maxWaveletSteps;

  // perform the wavelet transform
  if (opt.waveletSteps > 0) {
    float waveletMs = haar_2d(size, data.floatData, false, opt.waveletSteps,
                              opt.isWaveletTransposeStandard);
    printf("Wavelet transform (%d steps): %.2f ms\n", 
           opt.waveletSteps, waveletMs);
  }

  // save the intermediate data to a file before quantizing
  if (opt.saveBeforeQuantizingFilename != "") {
    const char *filename = opt.saveBeforeQuantizingFilename.c_str();
    if (!writeDataFile(filename, data.floatData, width, height)) {
                       
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
  opt.thresholdValue = thresh_cpu(count, data.floatData, opt.thresholdFraction,
                                  &nonzeroCount, &maxVal, &minVal,
                                  &sortedAbsData);
			       
  maxAbsVal = sortedAbsData[count-1];

  // NOTE: nonzeroData points to data in sortedAbsData, don't deallocate both
  float *nonzeroData = sortedAbsData + count - nonzeroCount;

  elapsed = NixTimer::time() - startTime;
  printf("Compute threshold = %g, min = %g, max = %g: %.2f ms\n",
	 opt.thresholdValue, minVal, maxVal, elapsed*1000);

  // std::vector<float> quantBinBoundaries;
  // std::vector<float> quantBinValues;

  quantizedData.initInts(width, height);

  // quantize the data
  startTime = NixTimer::time();
  switch (opt.quantizeAlgorithm) {

  case QUANT_ALG_UNIFORM:
    {
      QuantUniform qunif(opt.quantizeBits, opt.thresholdValue, maxAbsVal);
      QuantizationLooper<QuantUniform> qloop(&qunif, opt.quantizeBits);
      qloop.quantize(count, data.floatData, quantizedData.intData, true);
      printf("Quantization error: %g\n", qloop.getError());
      opt.maxAbsVal = maxAbsVal;
    }
    break;

  case QUANT_ALG_LOG:
    {
      QuantLog qlog(opt.quantizeBits, opt.thresholdValue, maxAbsVal);
      QuantizationLooper<QuantLog> qloop(&qlog, opt.quantizeBits);
      qloop.quantize(count, data.floatData, quantizedData.intData, true);
      printf("Quantization error: %g\n", qloop.getError());
      opt.maxAbsVal = maxAbsVal;
    }
    break;

  case QUANT_ALG_COUNT:
    {
      QuantCodebook qcb;
      qcb.initCountBins(count, data.floatData, opt.quantizeBits,
                        opt.thresholdValue);
      // qcb.printCodebook();
      opt.quantBinBoundaries = qcb.boundaries;
      opt.quantBinValues = qcb.codebook;
      QuantizationLooper<QuantCodebook> qloop(&qcb, opt.quantizeBits);
      qloop.quantize(count, data.floatData, quantizedData.intData, true);
      printf("Quantization error: %g\n", qloop.getError());
    }
    break;

  case QUANT_ALG_LLOYD: 
    {
      computeLloydQuantization(nonzeroData, nonzeroCount,
			       opt.quantizeBits,
			       opt.quantBinBoundaries, opt.quantBinValues);
      elapsed = NixTimer::time() - startTime;
      printf("Lloyd quantization %.2f ms\n", elapsed*1000);
      startTime = NixTimer::time();
      QuantCodebook qcb(opt.quantBinBoundaries, opt.quantBinValues);
      // qcb.printCodebook();
      QuantizationLooper<QuantCodebook> qloop(&qcb, opt.quantizeBits);
      qloop.quantize(count, data.floatData, quantizedData.intData, true);
      printf("Quantization error: %g\n", qloop.getError());
    }
    break;

  default:
    fprintf(stderr, "Quantization algorithm %d not found.\n",
            (int)opt.quantizeAlgorithm);
    return false;
  }

  elapsed = NixTimer::time() - startTime;
  printf("Quantize: %.2f ms\n", elapsed*1000);

  nonzeroData = NULL;
  delete[] sortedAbsData;
  sortedAbsData = NULL;

  // Try out huffman encoding
  // testHuffman(quantizedData, count, opt.quantizeBits);

  // write the quantized data to a file

  startTime = NixTimer::time();

  if (!writeQuantData(outputFile, quantizedData, opt))
    return false;
    
  elapsed = NixTimer::time() - startTime;
  printf("Write data file: %.2f ms\n", elapsed*1000);

  elapsed = NixTimer::time() - firstStartTime;
  printf("Total: %.2f ms\n", elapsed*1000);

  return true;
}

bool decompressFile(const char *inputFile, const char *outputFile,
		    Options &opt) {

  Data2d quantizedData, data;
  double firstStartTime, startTime, elapsed;

  firstStartTime = startTime = NixTimer::time();

  if (!readQuantData(inputFile, quantizedData, opt)) return false;

  const int width = quantizedData.width;
  const int height = quantizedData.height;

  elapsed = NixTimer::time() - startTime;
  printf("Read %dx%d data file: %.2f ms\n", width, height,
	 (NixTimer::time() - startTime) * 1000);

  assert(quantizedData.intData != NULL && width > 0 && height > 0);
  assert(width == height);
  int size = width;
  int count = width * height;
  data.initFloats(width, height);

  // de-quantize the data
  startTime = NixTimer::time();
  switch (opt.quantizeAlgorithm) {
  case QUANT_ALG_UNIFORM:
    {
      QuantUniform qunif(opt.quantizeBits, opt.thresholdValue, opt.maxAbsVal);
      QuantizationLooper<QuantUniform> qloop(&qunif, opt.quantizeBits);
      qloop.dequantize(count, quantizedData.intData, data.floatData);
    }
    break;

  case QUANT_ALG_LOG:
    {
      QuantLog qunif(opt.quantizeBits, opt.thresholdValue, opt.maxAbsVal);
      QuantizationLooper<QuantLog> qloop(&qunif, opt.quantizeBits);
      qloop.dequantize(count, quantizedData.intData, data.floatData);
    }
    break;

  case QUANT_ALG_COUNT:
  case QUANT_ALG_LLOYD:
    {
      QuantCodebook qcb;
      qcb.init(opt.quantBinBoundaries, opt.quantBinValues);
      QuantizationLooper<QuantCodebook> qloop(&qcb, opt.quantizeBits);
      qloop.dequantize(count, quantizedData.intData, data.floatData);
    }
    
    break;

  default:
    fprintf(stderr, "Quantization algorithm %d not found.\n",
            (int)opt.quantizeAlgorithm);
    return false;
  }

  printf("Dequantize: %.2f ms\n", (NixTimer::time() - startTime) * 1000);

  // perform inverse wavelet transform
  elapsed = haar_2d(size, data.floatData, true, opt.waveletSteps,
                    opt.isWaveletTransposeStandard);
  printf("Wavelet inverse transform: %.2f ms\n", elapsed);

  // write the reconstructed data
  startTime = NixTimer::time();
  if (!writeDataFile(outputFile, data.floatData, data.width, data.height))
    return false;
  printf("Write file: %.2f ms\n", (NixTimer::time() - startTime) * 1000);

  printf("Total: %.2f ms\n", (NixTimer::time() - firstStartTime) * 1000);

  return true;
}


/*
  Distribute N codebook entries like this:

                                      1
                                      |
       N/2-1        1        N/2-1    v
  ---negatives--+-------+--positives--x
                ^   0   ^             ^
                |       |             |
          -thresh       +thresh      max

  Except for the maximum positive value, which has its own codebook
  entry, negative and positive entries will be mirrored.

  For example, if N = 8, thresh=1, and max=10, one possible
  set of positive codebook entries is: 1, 3, 7
  All codebook entries:
   -7, -3, 1, 0, 1, 3, 7, 10
*/
void computeLloydQuantization
(const float *inputData, int count, 
 int bits,  // # of quantization bits
 std::vector<float> &quantBinBoundaries,
 std::vector<float> &quantBinValues) {
  
  int binCount = (1 << (bits - 1)) - 1;
  quantBinBoundaries.clear();
  quantBinValues.clear();
  float *binBoundaries = new float[binCount-1];
  float *binValues = new float[binCount];

  // inputData is sorted, use it to get minVal and maxVal
  float minVal = inputData[0];

  // Skip the last entry in inputData[] because it is often much larger than
  // the rest and skews the results
  float maxVal = inputData[count-2];

  assert(minVal > 0);
  assert(maxVal > minVal);

  /* use log quantization to create an initial codebook
     f(minVal) = 0, f(maxVal) = binCount
     f(x) = b*log(a*x)

       b*log(a*min) = 0  b*log(a*max) = binCount
       log(a*min) = 0    b = binCount / log(a*max)
       a*min = 1
       a = 1/min


     y = b*log(a*x)
     y/b = log(a*x)
     e^(y/b) = a*x
     e^(y/b) / a = x

       1/a = min, logScale = 1/b = log(max/min) / binCount

     min * e^(y*logScale) = x
  */

  // printf("min=%f, max=%f\n", minVal, maxVal);
  float logScale = logf(maxVal / minVal) / binCount;
  for (int i=0; i < binCount; i++) {
    binValues[i] = minVal * expf(i * logScale);
  }

  /*
  for (int i=0; i < binCount; i++) {
    printf("InitCB %d: %f\n", i, binValues[i]);
  }
  */

  // fine-tune the codebook and bin boundaries using Lloyd's algorithm.
  // This also applies the quantization to each value, writing the values
  // to quantizedData[].
  float dist, reldist;
  unsigned *quantizedData = new unsigned[count];
  lloyd(inputData, count-1, binValues, binCount, binBoundaries, dist,
        reldist, quantizedData);

  /*
  printf("LLoyd output:\n");
  for (int i=0; ; i++) {
    printf("codebook %d: %g\n", i, binValues[i]);
    if (i == binCount-1) break;
    printf("bin %d: %g\n", i, binBoundaries[i]);
  }
  */

  // sanity-check
  for (int i=0; i < binCount-1; i++) {
    if (binValues[i] > binValues[i+1]) {
      printf("ERROR: codebook[%d] > codebook[%d]  (%f > %f)\n",
             i, i+1, binValues[i], binValues[i+1]);
    }

    if (binBoundaries[i] < binValues[i] || binBoundaries[i] > binValues[i+1]) {
      printf("ERROR: partition[%d] (%.8g) should be between codebook[%d] (%.8g) and codebook[%d] (%.8g)\n",
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
  quantBinBoundaries.push_back(inputData[count-1]);
  quantBinValues.push_back(inputData[count-1]);

  /*
    // print all the input data and the quantized value for each
  for (int i=0; i < count; i++) {
    printf("%g\t%d\n", inputData[i], quantizedData[i]);
  }
  */
  
  delete[] quantizedData;
  delete[] binBoundaries;
  delete[] binValues;
}

  
void testHuffman(const int data[], int count, int bitCount) {
  double startTime = NixTimer::time();
  int valueCount = 1 << bitCount;
  Huffman huff(valueCount);
  
  for (int i=0; i < count; i++) {
    huff.increment(data[i]);
  }

  huff.computeHuffmanCoding();
  double elapsed = NixTimer::time() - startTime;
  printf("Huffman build table %.3f ms\n", elapsed*1000);

  // huff.printEncoding();
  
  startTime = NixTimer::time();
  FILE *f = fopen("huff.out", "wb");
  BitStreamWriter bitWriter(f);
  huff.encodeToStream(&bitWriter, data, count);
  bitWriter.flush();
  fclose(f);
  // printf("%llu bits written\n", (long long unsigned) bitWriter.size());
  size_t bitsWritten = bitWriter.size();
  int bytesWritten = (bitsWritten + 31) / 32 * 4;

  long long unsigned totalBits = 0;
  for (int i=0; i < valueCount; i++) {
    totalBits += huff.encodedLength(i) * huff.getCount(i);
  }

  printf("Huffman encoding: %d bytes, %.2f bits/pixel, "
	 "longest encoding = %d bits\n",
	 bytesWritten, (double)totalBits / count,
	 huff.getLongestEncodingLength());
  elapsed = NixTimer::time() - startTime;
  printf("Huffman write file %.3f ms\n", elapsed*1000);

  startTime = NixTimer::time();
  f = fopen("huff.out", "rb");
  BitStreamReader bitReader(f);
  int *values2 = new int[count];
  int readCount = huff.decodeFromStream(values2, count, &bitReader);
  fclose(f);
  elapsed = NixTimer::time() - startTime;
  printf("Huffman read file %.3f ms\n", elapsed*1000);

  if (count != readCount)
    printf("ERROR wrote %d encoded values, read %d\n", count, readCount);

  for (int i=0; i < count; i++) {
    if (values2[i] != data[i]) {
      printf("Error at %d: got %d, should be %d\n", i, values2[i], data[i]);
      break;
    }
  }
  delete[] values2;
}
