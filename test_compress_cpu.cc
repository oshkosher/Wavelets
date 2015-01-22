#include <google/protobuf/stubs/common.h>
#include <algorithm>
#include "test_compress_common.h"
#include "dwt_cpu.h"
#include "nixtimer.h"
#include "thresh_cpu.h"
#include "quant.h"
#include "lloyds.h"
#include "bit_stream.h"
#include "huffman.h"
#include "wavelet.h"
#include "cubelet_file.h"

using namespace scu_wavelet;

// global variable that disables status output
bool QUIET = false;

// returns true iff 'suffix' is a suffix of 's'
bool endsWith(const char *s, const char *suffix);
bool compressFile(const char *inputFile, const char *outputFile, Options &opt);
bool decompressFile(const char *inputFile, const char *outputFile,
		    Options &opt);

// make the cubelet data into floats, if it isn't already
void translateCubeDataToFloat(Cube *src, CubeFloat *dest);

// turn the data back into the original datatype, if it wasn't floats
void translateCubeFloatToData(CubeFloat *src, Cube *dest, bool verbose = false);


void computeLloydQuantization(const float *inputData, int count, int bits,
			      std::vector<float> &quantBinBoundaries,
			      std::vector<float> &quantBinValues);
void testHuffman(const int data[], int count, int bitCount);

// Read one cubelet from the input file into 'inputData' (the original data)
// and a copy in 'data' (padded and translated to floats).
bool readData(const char *filename,
              Cube *inputData,  // the original data, may not be floats
              CubeFloat *data,  // input data transformed into floats
              CubeletStreamReader *cubeletStream,
              WaveletCompressionParam *param);

// pad the data so each transform step has an even number of elements
void padDataSize(int3 &size, int3 steps);

// pad the actual data, repeating the last value in each axis as needed
void padData(CubeFloat &data, int3 originalSize);

void setWaveletSteps(int3 &steps, const int3 &size);

// Peform the wavelet transform
bool waveletTransform(CubeFloat &data, const WaveletCompressionParam &param,
                      bool isInverse, bool verbose);

// quantize the data
bool quantize(const CubeFloat &data, CubeInt &quantizedData,
              float maxAbsVal, WaveletCompressionParam &param,
              const float *nonzeroData, int nonzeroCount,
              float *quantErrorOut = NULL);

bool dequantize(const CubeInt &quantizedData, CubeFloat &data,
                const WaveletCompressionParam &param);

void computeErrorRates(const CubeInt *quantizedData,
                       const WaveletCompressionParam &param,
                       const Cube *inputData,
                       float *meanSquaredError,
                       float *peakSNR);

void quantizationExperiments(CubeFloat &data, Options &opt,
                             const float *sortedAbsData,
                             float maxAbsVal);

/*
//  test ErrorAccumulator
int main() {
  ErrorAccumulator err;
  err.add(10, 12);
  err.add(20, 17);
  err.add(30, 35);
  err.add(40, 50);
  err.setMaxPossible(255);
  printf("L1 %f, L2 %f, PSNR %f\n", err.getAverageError(),
         err.getMeanSquaredError(), err.getPeakSignalToNoiseRatio());
  return 0;
}
*/
  

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


// returns true iff 'suffix' is a suffix of 's'
bool endsWith(const char *s, const char *suffix) {
  int len = strlen(s), suffixLen = strlen(suffix);
  if (suffixLen > len) return false;
  return 0 == strcmp(s + len - suffixLen, suffix);
}


bool compressFile(const char *inputFile, const char *outputFile,
                  Options &opt) {

  WaveletCompressionParam param = opt.param;
  Cube inputData;
  CubeFloat data;
  CubeInt quantizedData;
  double firstStartTime, startTime, elapsed;
  CubeletStreamReader cubeletStream;

  firstStartTime = NixTimer::time();

  // read and pad the data
  if (!readData(inputFile, &inputData, &data, &cubeletStream, &param))
    return false;

  // perform the wavelet transformation
  if (!waveletTransform(data, param, false, opt.verbose)) return false;

  // save the intermediate data to a file before quantizing
  if (opt.saveBeforeQuantizingFilename != "") {
    const char *filename = opt.saveBeforeQuantizingFilename.c_str();
    CubeletStreamWriter writer;
    data.param = param;

    // if the data is 1xHxD, transpose it to HxDx1
    bool isTransposed = false;
    if (data.width() == 1) {
      data.transpose3dFwd();
      isTransposed = true;
    }
    if (!writer.open(filename) ||
        !writer.addCubelet(&data) ||
        !writer.close()) {
      fprintf(stderr, "Failed to write intermediate data file \"%s\".\n",
              filename);
    } else {
      if (!QUIET)
        printf("Write intermediate data file \"%s\"\n", filename);
    }
    if (isTransposed)
      data.transpose3dBack();
  }


  // find the threshold value by sorting
  // XXX quickselect will speed up this selection, but then we'll lose the
  // benefits of having the sorted data
  const int count = data.count();
  float maxVal, minVal, maxAbsVal, *sortedAbsData = NULL;
  int nonzeroCount;

  startTime = NixTimer::time();

  param.thresholdValue
    = thresh_cpu(&data, param.thresholdFraction, &nonzeroCount, &maxVal,
                 &minVal, &sortedAbsData);

  maxAbsVal = sortedAbsData[count-1];

  // nonzeroCount is the number of nonzero entries (those whose absolute
  // values are greater than thresholdValue) in sortedAbsData[].
  // Since sortedAbsData is sorted, those are all at the end of the array.
  // Set nonzeroData to point to the beginning of that nonzero data.
  // NOTE: nonzeroData points to data in sortedAbsData, don't deallocate both
  float *nonzeroData = sortedAbsData + count - nonzeroCount;

  elapsed = NixTimer::time() - startTime;
  if (!QUIET)
    printf("Compute threshold = %g, min = %g, max = %g: %.2f ms\n",
           param.thresholdValue, minVal, maxVal, elapsed*1000);

  // don't write a data file; just run some experiments testing different
  // quantization settings
  /*
  if (opt.runQuantizationExperiments) {
    quantizationExperiments(data, opt, sortedAbsData, maxAbsVal);
    return true;
  }
  */

  if (!quantize(data, quantizedData, maxAbsVal, param, nonzeroData,
                nonzeroCount)) return false;

  // compute error rates
  if (opt.doComputeError) {
    float meanSqErr, peakSNR;
    computeErrorRates(&quantizedData, param, &inputData, &meanSqErr, &peakSNR);
    printf("Mean squared error: %.3f, peak SNR: %.3f\n", meanSqErr, peakSNR);
  }

  if (opt.verbose) quantizedData.print("After quantization");
  
  // deallocate sortedAbsData
  delete[] sortedAbsData;
  nonzeroData = NULL;
  sortedAbsData = NULL;

  // testHuffman(quantizedData.intData, quantizedData.count(), opt.quantizeBits);

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

  return true;
}


bool decompressFile(const char *inputFile, const char *outputFile,
		    Options &opt) {

  CubeInt quantizedData;
  CubeFloat data;
  Cube outputData;
  double firstStartTime, startTime, elapsed;

  firstStartTime = startTime = NixTimer::time();

  // read the compressed file
  CubeletStreamReader cubeletStream;
  if (!cubeletStream.open(inputFile)) return false;
  if (!readQuantData(cubeletStream, &quantizedData)) return false;

  WaveletCompressionParam &param = quantizedData.param;

  elapsed = NixTimer::time() - startTime;
  if (!QUIET)
    printf("Read %dx%dx%d data file: %.2f ms\n",
           quantizedData.size.x, quantizedData.size.y, quantizedData.size.z,
           elapsed * 1000);

  data.size = quantizedData.size;
  data.param = param;
  data.allocate();

  if (opt.verbose) quantizedData.print("Before dequantize");
  
  // de-quantize the data
  dequantize(quantizedData, data, param);

  // perform inverse wavelet transform
  if (!waveletTransform(data, param, true, opt.verbose)) return false;

  // change the data back to the original datatype, if necessary
  outputData.size = param.originalSize;
  outputData.datatype = param.originalDatatype;
  translateCubeFloatToData(&data, &outputData, opt.verbose);

  // write the reconstructed data
  startTime = NixTimer::time();
  CubeletStreamWriter outStream;
  if (!outStream.open(outputFile) ||
      !outStream.addCubelet(&outputData) ||
      !outStream.close())
    return false;

  if (!QUIET) {
    printf("Write file: %.2f ms\n", (NixTimer::time() - startTime) * 1000);
    printf("Total: %.2f ms\n", (NixTimer::time() - firstStartTime) * 1000);
  }

  return true;
}


/**
   Read a cubelet from the input file.
     inputData will contain the original data
     data will contain a copy of the data, padded appropriately for the
       the number of transform steps in param->transformSteps
       (for n steps, pad up to a multiple of 2^n)
     How param is used:
       transformSteps is read to do the padding
       originalSize is set
       originalDatatype is set

   If the filename ends with ".cube":
     It is assumed to be a cubelet file, and is read with
     cubeletStream. cubeletStream will be opened if it isn't
     already. It will not be closed, so that more cubelets can be read.

   If the filename doesn't end with ".cube":
     It will be read with "readDataFile" from data_io.h.

   On error, output an error message and return false
*/
bool readData(const char *filename,
              Cube *inputData,  // the original data, may not be floats
              CubeFloat *data,  // input data transformed into floats
              CubeletStreamReader *cubeletStream,
              WaveletCompressionParam *param) {

  double startTime = NixTimer::time();

  bool isCubeletFile = endsWith(filename, ".cube");

  // get the size of the data without reading the data
  if (!isCubeletFile) {

    if (!readDataFile(filename, (float**)NULL,
                      &inputData->size.x, &inputData->size.y))
      return false;
    inputData->size.z = 1;
    inputData->datatype = WAVELET_DATA_FLOAT32;

  } else {

    // open the stream if it isn't already
    if (!cubeletStream->isOpen()) {
      if (!cubeletStream->open(filename)) return false;
    }

    // find the next uncompressed cubelet
    do {
      if (!cubeletStream->next(inputData)) {
        fprintf(stderr, "No uncompressed cubelets in the input data.\n");
        return false;
      }

    } while (inputData->isWaveletCompressed);

    // note: inputData.datatype might not be FLOAT32
  }

  inputData->allocate();
  
  // save the original size of the input data
  param->originalSize = inputData->size;
  param->originalDatatype = inputData->datatype;

  // adjust the number of wavelet steps in case the user requested too many
  setWaveletSteps(param->transformSteps, inputData->size);

  // Pad the size of the data (without touching the data itself, because
  // it hasn't been loaded yet) to align with the number of transformation
  // steps we'll be doing.
  int3 paddedSize = inputData->size;
  padDataSize(paddedSize, param->transformSteps);
  data->size = paddedSize;
  data->allocate();

  // read the actual data
  if (!isCubeletFile) {
    assert(inputData->datatype == WAVELET_DATA_FLOAT32);

    float *array;
    int width, height;
    if (!readDataFile(filename, &array, &width, &height)) return false;

    assert(inputData->size == int3(width, height, 1));
    memcpy(inputData->data_, array, sizeof(float)*width*height);
    data->copyFrom(array, inputData->size);
    delete[] array;

  } else { // read cubelet entry

    if (!cubeletStream->getCubeData(inputData)) return false;

    // translate the cube data if necessary
    translateCubeDataToFloat(inputData, data);

  }

  double elapsed = NixTimer::time() - startTime;
  if (!QUIET)
    printf("Read %dx%dx%d data file: %.2f ms\n",
           inputData->size.x, inputData->size.y, inputData->size.z,
           elapsed * 1000);

  padData(*data, inputData->size);

  return true;
}


class RowByteToFloat {
  CubeFloat *dest;
public:
  RowByteToFloat(CubeFloat *d) : dest(d) {}

  void visitRow(const unsigned char *readp, int length, int y, int z) {
    const unsigned char *end = readp + length;
    float *writep = dest->pointer(0, y, z);
    while (readp < end) {
      *writep++ = *readp++ * (1.0f / 255) - 0.5f;
    }
  }
};

class RowIntToFloat {
  CubeFloat *dest;
public:
  RowIntToFloat(CubeFloat *d) : dest(d) {}

  void visitRow(const int *readp, int length, int y, int z) {
    const int *end = readp + length;
    float *writep = dest->pointer(0, y, z);
    while (readp < end) {
      *writep++ = *readp++;
    }
  }
};

// make the cubelet data into floats, if it isn't already
void translateCubeDataToFloat(Cube *src, CubeFloat *dest) {

  // check that there are no insets
  assert(src->size == src->totalSize);

  // copy source cubelet id
  dest->parentOffset = src->parentOffset;

  if (src->datatype == WAVELET_DATA_FLOAT32) {
    dest->copyFrom(*(const CubeFloat*)src);
    return;
  }

  // check that dest's memory is already allocated
  assert(dest->data_);

  // int count = src->count();
  // float *writep = dest->pointer(0,0,0);
  
  // map 0..255 to -.5 .. .5
  if (src->datatype == WAVELET_DATA_UINT8) {
    CubeByte *s = (CubeByte*) src;
    RowByteToFloat rowVisitor(dest);
    s->visitRows(rowVisitor);
  }

  // without any other range information, just translate directly to floats
  else if (src->datatype == WAVELET_DATA_INT32) {
    CubeInt *s = (CubeInt*) src;
    RowIntToFloat rowVisitor(dest);
    s->visitRows(rowVisitor);
  }
}


// turn the data back into the original datatype, if it wasn't floats
// dest->size has been set, and it may be smaller than src
void translateCubeFloatToData(CubeFloat *src, Cube *dest, bool verbose) {

  // check that there are no insets
  assert(src->size == src->totalSize);

  // copy source size
  dest->totalSize = src->totalSize;
  dest->parentOffset = src->parentOffset;

  if (dest->datatype == src->datatype) {
    dest->data_ = src->data_;
    src->data_ = NULL;
    dest->ownsData = src->ownsData;
    return;
  }

  dest->allocate();

  int count = src->count();
  float *readp = src->pointer(0,0,0), *endp = readp + count;
  
  // map -.5 .. .5 to 0..255
  if (dest->datatype == WAVELET_DATA_UINT8) {
    CubeByte *d = (CubeByte*) dest;
    unsigned char *writep = d->pointer(0,0,0);

    while (readp < endp) {
      int x = (*readp++ + 0.5f) * 255;
      if (x < 0) {
        x = 0;
      } else if (x > 255) {
        x = 255;
      }
      *writep++ = x;
    }
    if (verbose) d->print("After restoring data type");
  }

  // without any other range information, just translate directly to ints
  else if (dest->datatype == WAVELET_DATA_INT32) {
    CubeInt *d = (CubeInt*) dest;
    int *writep = d->pointer(0,0,0);

    // map 0..255 to -.5 .. .5
    while (readp < endp) {
      *writep++ = *readp++;
    }

    if (verbose) d->print("After restoring data type");
  }

  src->deallocate();

}


// pad the data so each transform step has an even number of elements
void padDataSize(int3 &size, int3 steps) {
  size.x = dwt_padded_length(size.x, steps.x, false);
  size.y = dwt_padded_length(size.y, steps.y, false);
  size.z = dwt_padded_length(size.z, steps.z, false);
}


// REPEAT-pad data
void padData(CubeFloat &data, int3 originalSize) {
  
  if (data.size == originalSize) return;
  assert(data.size >= originalSize);

  double startTime = NixTimer::time();

  // for each existing layer of data, pad the end of each row and column
  for (int z = 0; z < originalSize.z; z++) {

    // pad each row by replicating the last value in each row
    for (int y = 0; y < originalSize.y; y++) {
      float value = data.get(originalSize.x-1, y, z);
      float *p = data.pointer(originalSize.x, y, z);
      float *e = p + data.width() - originalSize.x;
      while (p < e) *p++ = value;
    }

    // pad each column by replicating the last row
    float *rowp = data.pointer(0, originalSize.y-1, z);
    for (int y = originalSize.y; y < data.height(); y++)
      memcpy(data.pointer(0, y, z), rowp, data.width() * sizeof(float));
  }

  // for the new layers of data, copy the last layer of data
  float *layerp = data.pointer(0, 0, originalSize.z-1);
  int layerSize = data.width() * data.height() * sizeof(float);
  for (int z=originalSize.z; z < data.depth(); z++)
    memcpy(data.pointer(0, 0, z), layerp, layerSize);

  double elapsed = NixTimer::time() - startTime;
  if (!QUIET)
    printf("Pad data from %dx%dx%d to %dx%dx%d: %.2f ms\n",
           originalSize.x, originalSize.y, originalSize.z, 
           data.width(), data.height(), data.depth(),
           elapsed*1000);

}


void setWaveletSteps(int3 &steps, const int3 &size) {
  int maxSteps = dwtMaximumSteps(size.x);
  if (steps.x < 0 || steps.x > maxSteps)
    steps.x = maxSteps;

  maxSteps = dwtMaximumSteps(size.y);
  if (steps.y < 0 || steps.y > maxSteps)
    steps.y = maxSteps;

  maxSteps = dwtMaximumSteps(size.z);
  if (steps.z < 0 || steps.z > maxSteps)
    steps.z = maxSteps;
}


// Peform the wavelet transform
bool waveletTransform(CubeFloat &data, const WaveletCompressionParam &param,
                      bool isInverse, bool verbose) {

  // do nothing if 0 steps in each direction
  if (param.transformSteps == int3(0,0,0)) return true;
  
  double startTime = NixTimer::time();

  if (verbose) data.print("Before wavelet transform");

  if (param.waveletAlg == WAVELET_CDF97) {
    cdf97_3d(&data, param.transformSteps, isInverse, 
             param.isWaveletTransposeStandard);
  }

  else if (param.waveletAlg == WAVELET_HAAR) {
    haar_3d(&data, param.transformSteps, isInverse, 
            param.isWaveletTransposeStandard);
  }

  else {
    fprintf(stderr, "Unknown wavelet: %s\n",
            waveletAlgToName(param.waveletAlg));
    return false;
  }

  if (verbose) data.print("After wavelet transform");

  double elapsedMs = (NixTimer::time() - startTime) * 1000;

  if (!QUIET)
    printf("%s%s wavelet transform (%d,%d,%d steps): %.2f ms\n", 
           waveletAlgToName(param.waveletAlg),
           isInverse ? " inverse" : "",
           param.transformSteps.x, param.transformSteps.y,
           param.transformSteps.z, elapsedMs);

  return true;
}


bool quantize(const CubeFloat &data, CubeInt &quantizedData,
              float maxAbsVal, WaveletCompressionParam &param,
              const float *nonzeroData, int nonzeroCount,
              float *quantErrorOut) {

  float quantErr = -1;

  // for now, to make the processing easier, don't accept padded data
  assert(data.size == data.totalSize);
  
  // XXX remove this when all the quantization algorithms switch from
  // quantizeBits to binCount
  int quantizeBits = ceilLog2(param.binCount);

  // allocate memory for the quantized data if it isn't already
  quantizedData.size = quantizedData.totalSize = data.size;
  quantizedData.inset = int3(0,0,0);
  quantizedData.allocate();

  const int count = data.count();
  const float *inputData = data.pointer(0,0,0);
  int *outputData = quantizedData.pointer(0,0,0);

  double startTime = NixTimer::time();

  switch (param.quantAlg) {
    
  case QUANT_ALG_UNIFORM:
    {  // use a new code block so we can allocate new variables
      QuantUniform qunif(quantizeBits, param.thresholdValue, maxAbsVal);
      QuantizationLooper<QuantUniform> qloop(&qunif, quantizeBits);
      qloop.quantize(count, inputData, outputData, quantErrorOut != NULL);
      quantErr = qloop.getError();
      param.maxValue = maxAbsVal;
    }
    break;

  case QUANT_ALG_LOG:
    {
      QuantLog qlog(quantizeBits, param.thresholdValue, maxAbsVal);
      QuantizationLooper<QuantLog> qloop(&qlog, quantizeBits);
      qloop.quantize(count, inputData, outputData, quantErrorOut != NULL);
      quantErr = qloop.getError();
      param.maxValue = maxAbsVal;
    }
    break;

  case QUANT_ALG_LLOYD: 
    {
      computeLloydQuantization(nonzeroData, nonzeroCount, quantizeBits,
			       param.binBoundaries, param.binValues);
      double elapsed = NixTimer::time() - startTime;
      if (!QUIET)
        printf("Lloyd quantization %.2f ms\n", elapsed*1000);
      startTime = NixTimer::time();
      QuantCodebook qcb(param.binBoundaries, param.binValues);
      // qcb.printCodebook();
      QuantizationLooper<QuantCodebook> qloop(&qcb, quantizeBits);
      qloop.quantize(count, inputData, outputData, quantErrorOut != NULL);
      quantErr = qloop.getError();
    }
    break;

  default:
    fprintf(stderr, "Quantization algorithm %d not found.\n",
            (int)param.quantAlg);
    return false;
  }

  // if the caller asked for the error rate to be computed print it and
  // return it
  if (quantErrorOut && quantErr >= 0) {
    if (!QUIET)
      printf("Quantization error: %g\n", quantErr);
    if (quantErrorOut) *quantErrorOut = quantErr;
  }

  if (!QUIET) {
    double elapsed = NixTimer::time() - startTime;
    printf("Quantize %s: %.2f ms\n", quantAlgId2Name(param.quantAlg),
           elapsed*1000);
  }

  return true;
}


bool dequantize(const CubeInt &quantizedData, CubeFloat &data,
                const WaveletCompressionParam &param) {

  double startTime = NixTimer::time();
  
  const int count = quantizedData.count();
  const int *inputData = quantizedData.pointer(0,0,0);
  float *outputData = data.pointer(0,0,0);

  // XXX remove this when all the quantization algorithms switch from
  // quantizeBits to binCount
  int quantizeBits = ceilLog2(param.binCount);

  switch (param.quantAlg) {
  case QUANT_ALG_UNIFORM:
    {
      QuantUniform qunif(quantizeBits, param.thresholdValue, param.maxValue);
      QuantizationLooper<QuantUniform> qloop(&qunif, quantizeBits);
      qloop.dequantize(count, inputData, outputData);
    }
    break;

  case QUANT_ALG_LOG:
    {
      QuantLog qunif(quantizeBits, param.thresholdValue, param.maxValue);
      QuantizationLooper<QuantLog> qloop(&qunif, quantizeBits);
      qloop.dequantize(count, inputData, outputData);
    }
    break;

  case QUANT_ALG_LLOYD:
    {
      QuantCodebook qcb;
      qcb.init(param.binBoundaries, param.binValues);
      QuantizationLooper<QuantCodebook> qloop(&qcb, quantizeBits);
      qloop.dequantize(count, inputData, outputData);
    }
    break;

  default:
    fprintf(stderr, "Quantization algorithm %d not found.\n",
            (int)param.quantAlg);
    return false;
  }

  if (!QUIET)
    printf("Dequantize %s: %.2f ms\n", quantAlgId2Name(param.quantAlg),
           (NixTimer::time() - startTime) * 1000);

  return true;
}


/**
   Compare the compressed data with the original data.
   Given quantized data, dequantize it and apply the inverse wavelet
   transform.
*/
void computeErrorRates(const CubeInt *quantizedData,
                       const WaveletCompressionParam &param,
                       const Cube *inputData,
                       float *meanSquaredError,
                       float *peakSNR) {

  if (inputData->datatype != WAVELET_DATA_UINT8) {
    fprintf(stderr, "Current implementation can only compute error rates if input data is unsigned bytes (0..255)\n");
    return;
  }

  // dequantize quantizedData into restoredData
  CubeFloat restoredData;
  restoredData.size = quantizedData->size;
  restoredData.allocate();
  if (!dequantize(*quantizedData, restoredData, param)) return;

  // perform the inverse wavelet transform on restoredData
  if (!waveletTransform(restoredData, param, true, false)) return;

  const int width = inputData->width();
  const CubeByte *original = (const CubeByte *) inputData;

  ErrorAccumulator errAccum;
  errAccum.setMaxPossible(255);  // largest unsigned char value

  // for each row of data in the original data, compare pixels
  for (int z=0; z < original->depth(); z++) {
    for (int y=0; y < original->height(); y++) {
      const unsigned char *originalRow = original->pointer(0, y, z);
      const float *restoredRow = restoredData.pointer(0, y, z);
      
      for (int x=0; x < width; x++) {
        float f = (restoredRow[x] + 0.5f) * 255;
        if (f < 0) {
          f = 0;
        } else if (f > 255) {
          f = 255;
        }
        unsigned char pixelValue = (unsigned char) f;

        errAccum.add(originalRow[x], pixelValue);
      }
    }
  }

  *meanSquaredError = errAccum.getMeanSquaredError();
  *peakSNR = errAccum.getPeakSignalToNoiseRatio();
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
  // Skip the last entry in inputData[] because it is often much larger than
  // the rest and skews the results
  float maxVal = inputData[count-2];
  assert(maxVal > 0);

  float minVal = inputData[0];
  if (minVal <= 0) {
    const float *nonzeroIdx = std::upper_bound(inputData, inputData+count, 0);
    minVal = *nonzeroIdx;
  }

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
      fprintf(stderr, "ERROR: codebook[%d] > codebook[%d]  (%f > %f)\n",
              i, i+1, binValues[i], binValues[i+1]);
    }

    if (binBoundaries[i] < binValues[i] || binBoundaries[i] > binValues[i+1]) {
      fprintf(stderr, "ERROR: partition[%d] (%.8g) should be between codebook[%d] (%.8g) and codebook[%d] (%.8g)\n",
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


/*
void quantizationExperiments(const Data2d &data, Options &opt,
                             const float *sortedAbsData,
                             float maxAbsVal) {
  
  Data2d quantizedData;
  Data3d data3d;
  Data3dInt quantizedData3d;

  // allocate memory once for the quantized data
  quantizedData.initInts(data.width, data.height);
  const int count = data.width * data.height;

  float quantErr;
  int fileSize;

  // print header row
  printf("Bits\t" "ThresholdPct\t" "Log10MeanSqErr\t" "CompressionRatio\n");

  // try different values for quantizeBits and thresholdFraction

  for (opt.quantizeBits = 4; opt.quantizeBits <= 12; opt.quantizeBits++) {

    for (float thresholdFraction = 0.15f; thresholdFraction < 0.95f;
         thresholdFraction += 0.05f) {
      
      int threshIdx = (int)(thresholdFraction * count);
      opt.thresholdValue = sortedAbsData[threshIdx];

      int nonzeroCount = count - threshIdx - 1;
      const float *nonzeroData = sortedAbsData + count - nonzeroCount;

      // quantize, saving the error rate
      if (!quantize(data, quantizedData, data3d, quantizedData3d,
                    maxAbsVal, opt,
                    nonzeroData, nonzeroCount, &quantErr)) break;

      // write to a dummy file, saving the file size
      if (!writeQuantData("test_compress_cpu.experiment.wz",
                          quantizedData, opt, &fileSize)) break;

      // check quantError just in case, to avoid calling log(0)
      if (quantErr <= 0) {
        quantErr = 0;
      } else {
        quantErr = log(quantErr) / log(10);
      }

      double compressionRatio = 100.0 * (count - fileSize) / count;

      // write a row of data
      printf("%d\t%.1f\t%.7g\t%f\n", opt.quantizeBits, thresholdFraction*100,
             quantErr, compressionRatio);
      fflush(stdout);
    }
  }
}
*/
