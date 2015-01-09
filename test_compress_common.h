#ifndef __TEST_COMPRESS_COMMON_H__
#define __TEST_COMPRESS_COMMON_H__

#include <vector>
#include <string>
#include "dwt_cpu.h"
#include "data_io.h"
#include "bit_stream.h"
#include "wavelet_compress.pb.h"

#define FILE_ID_STRING "SCU wavelet 1.0\n"

typedef enum {
  QUANT_ALG_UNKNOWN = -1,
  QUANT_ALG_UNIFORM = 1,
  QUANT_ALG_LOG = 2,
  QUANT_ALG_COUNT = 3,
  QUANT_ALG_LLOYD = 4
} QuantizeAlgorithm;


#define DEFAULT_WAVELET_STEPS -1
#define DEFAULT_THRESHOLD_FRACTION .5
#define DEFAULT_QUANTIZE_BITS 8
#define DEFAULT_QUANTIZE_ALGORITHM QUANT_ALG_LLOYD
#define DEFAULT_WAVELET_TRANSPOSE_STANDARD false
#define DEFAULT_PADDING_METHOD REFLECT

/** This holds the the parameters as set by the user on the command line,
    and the ones saved in the data file. */
struct Options {

  // if false, then decompress
  bool doCompress;

  // suppress output
  bool quiet;

  // # of wavelet steps. Negative value does the maximum.
  int waveletSteps;

  // if true, then do all wavelet steps in the x direction before doing all
  //          steps in the y direction
  // if false, do one x step, then one y step, and repeat
  bool isWaveletTransposeStandard;

  // fraction of the data that is rounded down to zero before quantizing
  float thresholdFraction;

  // # of bits into which the data is quantized
  int quantizeBits;

  // the algorithms used to quantize
  QuantizeAlgorithm quantizeAlgorithm;

  // After computing the bit encoding of each value, print the encodings
  bool printHuffmanEncoding;

  // -bq option : if not "", save a copy data in this file before quantizing
  std::string saveBeforeQuantizingFilename;


  // the following data is just for the save file; it does not correspond
  // to command line arguments

  /*
  int width, height;  // size of the 2-d image
  float *floatData;   // pre-quantized data is stored here
  int *intData;       // post-quantized data is stored here
  */

  float thresholdValue; // value computed via thresholdFraction

  // maximum absolute value in the data; used in uniform and log quantization
  float maxAbsVal;  

  // boundaries between codebook entries, use this to quantize data
  std::vector<float> quantBinBoundaries;

  // codebook - values assigned to each quantized value, use this to dequantize
  std::vector<float> quantBinValues;

  // table use to decode Huffman-encoded bits. See huffman.h (HuffmanDecoder)
  // for an explanation
  std::vector<int> huffmanDecodeTable;

  bool runQuantizationExperiments;

  Options() {
    init();
  }

  void init() {
    doCompress = true;
    quiet = false;
    waveletSteps = DEFAULT_WAVELET_STEPS;
    isWaveletTransposeStandard = DEFAULT_WAVELET_TRANSPOSE_STANDARD;
    thresholdFraction = DEFAULT_THRESHOLD_FRACTION;
    quantizeBits = DEFAULT_QUANTIZE_BITS;
    quantizeAlgorithm = DEFAULT_QUANTIZE_ALGORITHM;
    printHuffmanEncoding = false;
    saveBeforeQuantizingFilename = "";

    /*
    width = height = -1;
    floatData = NULL;
    intData = NULL;
    */

    thresholdValue = 0;
    maxAbsVal = 0;
    quantBinBoundaries.clear();
    quantBinValues.clear();
    huffmanDecodeTable.clear();

    runQuantizationExperiments = false;
  }
};


// Store 2d data
struct Data2d {
  int width, height;
  float *floatData;
  int *intData;

  Data2d() {
    width = height = -1;
    floatData = NULL;
    intData = NULL;
  }

  ~Data2d() {
    if (floatData) delete[] floatData;
    if (intData) delete[] intData;
  }

  void initInts(int width_, int height_) {
    assert(intData == NULL);
    width = width_;
    height = height_;
    floatData = NULL;
    intData = new int[width*height];
  }

  void initFloats(int width_, int height_) {
    assert(floatData == NULL);
    width = width_;
    height = height_;
    floatData = new float[width*height];
    intData = NULL;
  }

  int count() const {return width*height;}
};



QuantizeAlgorithm quantAlgName2Id(const char *name);
const char *quantAlgId2Name(QuantizeAlgorithm id);
WaveletCompressedImage_QuantizationAlgorithm quantAlgId2ProtoId
  (QuantizeAlgorithm id);
QuantizeAlgorithm quantProtoId2AlgId
  (WaveletCompressedImage_QuantizationAlgorithm protoId);

void printHelp();
bool parseOptions(int argc, char **argv, Options &opt, int &nextArg);

// Write fileData to a file
// If fileSizeBytes is not NULL, store the size of the output file in it.
bool writeQuantData(const char *filename, const Data2d &data, Options &opt,
                    int *fileSizeBytes = NULL);

// Read FileData from a file
bool readQuantData(const char *filename, Data2d &data, Options &opt);


#endif // __TEST_COMPRESS_COMMON_H__
