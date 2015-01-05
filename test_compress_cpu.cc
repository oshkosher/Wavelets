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
      QuantizationLooper<QuantUniform> qloop(&qunif, opt.quantizeBits);
      quantizedData = new int[size*size];
      qloop.quantize(size*size, data, quantizedData, true);
      printf("Quantization error: %g\n", qloop.getError());
    }
    break;

  case QUANT_ALG_LOG:
    // quant_log_cpu(size, data, opt.quantizeBits, threshold, maxVal);
    {
      QuantLog qlog(opt.quantizeBits, threshold, maxAbsVal);
      QuantizationLooper<QuantLog> qloop(&qlog, opt.quantizeBits);
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
      // qcb.printCodebook();
      quantBinBoundaries = qcb.boundaries;
      quantBinValues = qcb.codebook;
      QuantizationLooper<QuantCodebook> qloop(&qcb, opt.quantizeBits);
      quantizedData = new int[count];
      qloop.quantize(count, data, quantizedData, true);
      printf("Quantization error: %g\n", qloop.getError());
    }
    break;

  case QUANT_ALG_LLOYD: 
    {
      computeLloydQuantization(nonzeroData, nonzeroCount,
			       opt.quantizeBits,
			       quantBinBoundaries, quantBinValues);
      elapsed = NixTimer::time() - startTime;
      printf("Lloyd quantization %.2f ms\n", elapsed*1000);
      startTime = NixTimer::time();
      QuantCodebook qcb(quantBinBoundaries, quantBinValues);
      // qcb.printCodebook();
      QuantizationLooper<QuantCodebook> qloop(&qcb, opt.quantizeBits);
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

  nonzeroData = NULL;
  delete[] sortedAbsData;
  sortedAbsData = NULL;

  // Try out huffman encoding
  // testHuffman(quantizedData, count, opt.quantizeBits);

  // write the quantized data to a file
  FileData fileData(opt, NULL, quantizedData, size, size);
  fileData.threshold = threshold;
  if (opt.quantizeAlgorithm == QUANT_ALG_UNIFORM ||
      opt.quantizeAlgorithm == QUANT_ALG_LOG) {
    fileData.quantMaxVal = maxVal;
  } else {
    fileData.quantBinBoundaries = quantBinBoundaries;
    fileData.quantBinValues = quantBinValues;
  }

  startTime = NixTimer::time();

  if (!writeQuantData(outputFile, fileData, opt.printHuffmanEncoding))
    return false;
    
  elapsed = NixTimer::time() - startTime;
  printf("Write data file: %.2f ms\n", elapsed*1000);

  elapsed = NixTimer::time() - firstStartTime;
  printf("Total: %.2f ms\n", elapsed*1000);

  if (quantizedData) delete[] quantizedData;

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

  if (!readQuantData(inputFile, f)) return false;

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
      QuantizationLooper<QuantUniform> qloop(&qunif, opt.quantizeBits);
      qloop.dequantize(count, inputData, data);
    }
    break;

  case QUANT_ALG_LOG:
    // dquant_log_cpu(size, data, f.quantizeBits, f.threshold, f.quantMaxVal);
    {
      QuantLog qunif(f.quantizeBits, f.threshold, f.quantMaxVal);
      QuantizationLooper<QuantLog> qloop(&qunif, opt.quantizeBits);
      qloop.dequantize(count, inputData, data);
    }
    break;

  case QUANT_ALG_COUNT:
  case QUANT_ALG_LLOYD:
    // dequant_codebook_array(f.quantBinValues, size*size, data);
    {
      QuantCodebook qcb;
      qcb.init(f.quantBinBoundaries, f.quantBinValues);
      QuantizationLooper<QuantCodebook> qloop(&qcb, opt.quantizeBits);
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
