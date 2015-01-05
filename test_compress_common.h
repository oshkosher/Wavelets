#ifndef __TEST_COMPRESS_COMMON_H__
#define __TEST_COMPRESS_COMMON_H__

#include <vector>
#include <string>
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


/** This holds the the parameters as set by the user on the command line. */
struct Options {
  bool doCompress;  // if false, then decompress
  int waveletSteps;
  bool isWaveletTransposeStandard;
  float thresholdFraction;
  int quantizeBits;
  bool printHuffmanEncoding; // print the bit encoding of each value
  std::string saveBeforeQuantizingFilename;  // -bq option
  QuantizeAlgorithm quantizeAlgorithm;
};

/** This holds all the data that will be written to or read from the data file. 
*/
struct FileData {
  float *floatData;
  int *intData;
  int width, height;

  int waveletSteps;
  int quantizeBits;
  float threshold;  // threshold value, not the proportion
  float quantMaxVal;  // maximum absolute value in the data;
                      // used in uniform and log quantization; 
  QuantizeAlgorithm quantizeAlgorithm;
  std::vector<float> quantBinBoundaries;
  std::vector<float> quantBinValues;

  std::vector<int> huffmanDecodeTable;

  // default constructor - invalid values
  FileData() : floatData(NULL), intData(NULL), width(-1), height(-1), 
    waveletSteps(0), quantizeBits(0), quantizeAlgorithm(QUANT_ALG_UNKNOWN) {}

  // alternative constructor - copy the options that make sense to
  // copy from 'Options'
  FileData(const Options &opt, float *floatData_=NULL, int *intData_=NULL,
           int width_=-1, int height_=-1) {
    floatData = floatData_;
    intData = intData_;
    width = width_;
    height = height_;
    waveletSteps = opt.waveletSteps;
    quantizeBits = opt.quantizeBits;
    quantizeAlgorithm = opt.quantizeAlgorithm;
    threshold = 0;
    quantMaxVal = 0;
  }
};


QuantizeAlgorithm quantAlgName2Id(const char *name);
const char *quantAlgId2Name(QuantizeAlgorithm id);
WaveletCompressedImage_QuantizationAlgorithm quantAlgId2ProtoId
  (QuantizeAlgorithm id);
QuantizeAlgorithm quantProtoId2AlgId
  (WaveletCompressedImage_QuantizationAlgorithm protoId);

class WriteRLEPairsToBitStream {
  BitStreamWriter *out;

  // # of bits in the value, not including sign bit
  int valueBits;

  // # of bits in the length
  int rlBits;
  
public:
  WriteRLEPairsToBitStream(BitStreamWriter *out_, int valueBits_, int rlBits_)
    : out(out_), valueBits(valueBits_), rlBits(rlBits_) {}
  
  void data(int value, int length) {
    // printf("%d * %d\n", length, value);

    // write sign bit
    out->write( (value<0) ? 1 : 0, 1);

    // write value
    out->write(abs(value), valueBits);
    
    // write length
    out->write(length, rlBits);
  }

  void end() {
    out->flush();
  }
};

void printHelp();
bool parseOptions(int argc, char **argv, Options &opt, int &nextArg);

// Write fileData to a file
bool writeQuantData(const char *filename, FileData &fileData,
                    bool printEncoding = false);

// Read FileData from a file
bool readQuantData(const char *filename, FileData &fileData);


#endif // __TEST_COMPRESS_COMMON_H__
