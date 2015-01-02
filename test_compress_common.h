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
  std::string saveBeforeQuantizingFilename;  // -bq option
  QuantizeAlgorithm quantizeAlgorithm;
};

/** This holds all the data that will be written to or read from the data
    file. We will be trying out a few different ways or writing the data
    file. All the routines will be given this same structure. */
struct FileData {
  float *data;
  int *intData;
  int width, height;

  int waveletSteps;
  int quantizeBits;
  float threshold;  // threshold value, not the proportion
  float quantMaxVal;
  QuantizeAlgorithm quantizeAlgorithm;
  std::vector<float> quantBinBoundaries;
  std::vector<float> quantBinValues;

  // default constructor - invalid values
  FileData() : data(NULL), intData(NULL), width(-1), height(-1), 
    waveletSteps(0), quantizeBits(0), quantizeAlgorithm(QUANT_ALG_UNKNOWN) {}

  // alternative constructor - copy the options that make sense to
  // copy from 'Options'
  FileData(const Options &opt, float *data_=NULL, int *intData_=NULL,
           int width_=-1, int height_=-1) {
    data = data_;
    intData = intData_;
    width = width_;
    height = height_;
    waveletSteps = opt.waveletSteps;
    quantizeBits = opt.quantizeBits;
    quantizeAlgorithm = opt.quantizeAlgorithm;
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
bool writeQuantDataSimple(const char *filename, FileData &fileData);
bool readQuantDataSimple(const char *filename, FileData &fileData);
bool writeQuantDataParamStrings(const char *filename, FileData &fileData);
bool readQuantDataParamStrings(const char *filename, FileData &fileData);
bool writeQuantDataProtoBuf(const char *filename, FileData &fileData);
bool readQuantDataProtoBuf(const char *filename, FileData &fileData);
bool readQuantData(const char *filename, FILE *inf, FileData &fileData);
bool writeQuantData(const char *filename, FILE *outf, FileData &fileData);


#endif // __TEST_COMPRESS_COMMON_H__
