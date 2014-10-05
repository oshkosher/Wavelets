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


#define DEFAULT_WAVELET_STEPS 3
#define DEFAULT_THRESHOLD_FRACTION .1
#define DEFAULT_QUANTIZE_BITS 8
#define DEFAULT_QUANTIZE_ALGORITHM QUANT_ALG_UNIFORM


/** This holds the the parameters as set by the user on the command line. */
struct Options {
  // alternative is to decompress
  bool doCompress;

  int waveletSteps;
  float thresholdFraction;
  int quantizeBits;
  QuantizeAlgorithm quantizeAlgorithm;
};

/** This holds all the data that will be written to or read from the data
    file. We will be trying out a few different ways or writing the data
    file. All the routines will be given this same structure. */
struct FileData {
  float *data;
  int width, height;

  int waveletSteps;
  int quantizeBits;
  float threshold;  // threshold value, not the proportion
  float quantMaxVal;
  QuantizeAlgorithm quantizeAlgorithm;
  std::vector<float> quantBinBoundaries;
  std::vector<float> quantBinValues;

  // default constructor - invalid values
  FileData() : data(NULL), width(-1), height(-1), waveletSteps(0),
		  quantizeBits(0), quantizeAlgorithm(QUANT_ALG_UNKNOWN) {}

  // alternative constructor - copy the options that make sense to
  // copy from 'Options'
  FileData(const Options &opt, float *data_=NULL, int width_=-1,
	      int height_=-1) {
    data = data_;
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
  
public:
  WriteRLEPairsToBitStream(BitStreamWriter *out_, int valueBits_)
    : out(out_), valueBits(valueBits_) {}
  
  // the length will always be 1..255
  void data(int value, int length) {
    // write sign bit
    out->write( (value<0) ? 1 : 0, 1);

    // write value
    out->write(abs(value), valueBits);
    
    // write length
    out->write(length, 8);
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
