/**
   Test the full path of tools:
     Read data file
     Perform wavelet transform
     Apply cutoff threshold
     Quantize
     Run-length encode
     Save to file

  And the inverse.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include "nixtimer.h"
#include "data_io.h"
#include "dwt_cpu.h"
#include "bit_stream.h"
#include "thresh_cpu.h"
#include "rle.h"
#include "quant_log_cpu.h"
#include "quant_unif_cpu.h"
#include "quant_count.h"
#include "dquant_log_cpu.h"
#include "dquant_unif_cpu.h"
#include "param_string.h"

using namespace std;

typedef enum {
  QUANT_ALG_UNKNOWN = -1,
  QUANT_ALG_UNIFORM = 1,
  QUANT_ALG_LOG = 2,
  QUANT_ALG_COUNT = 3,  // XXX not fully implemented
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
  vector<float> quantBinBoundaries;
  vector<float> quantBinValues;

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

void printHelp();
bool parseOptions(int argc, char **argv, Options &opt, int &nextArg);
bool compressFile(const char *inputFile, const char *outputFile, Options &opt);
bool decompressFile(const char *inputFile, const char *outputFile,
		    Options &opt);
bool writeQuantDataSimple(const char *filename, FileData &fileData);
bool readQuantDataSimple(const char *filename, FileData &fileData);
bool writeQuantDataParamStrings(const char *filename, FileData &fileData);
bool readQuantDataParamStrings(const char *filename, FileData &fileData);
bool readQuantData(const char *filename, FILE *inf, FileData &fileData);
bool writeQuantData(const char *filename, FILE *outf, FileData &fileData);

QuantizeAlgorithm quantAlgName2Id(const char *name);
const char *quantAlgId2Name(QuantizeAlgorithm id);

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


void printHelp() {
  printf("\n"
         "  test_compress [opts] <input file> <output file>\n"
         "  Options:\n"
	 "    -d : decompress (default is to compress)\n"
         "    -steps <n> : # of wavelet transform steps (default = %d)\n"
         "    -thresh <n> : proportion of data to be eliminated in threshold cutoff\n"
         "                  Must be between [0..1]. (default = %.3f)\n"
         "    -qbits <n> : # of bits into which data is quantized. Must be between \n"
         "                 1 and 32.  (default=%d)\n"
         "    -qalg <alorithm> : quantization algorithm: uniform, log, count, or lloyd\n"
         "                       (default = %s)\n"
         "\n",
         DEFAULT_WAVELET_STEPS,
         DEFAULT_THRESHOLD_FRACTION,
         DEFAULT_QUANTIZE_BITS,
         quantAlgId2Name(DEFAULT_QUANTIZE_ALGORITHM)
         );
  exit(1);
}


bool parseOptions(int argc, char **argv, Options &opt, int &nextArg) {
  opt.doCompress = true;
  opt.waveletSteps = DEFAULT_WAVELET_STEPS;
  opt.thresholdFraction = DEFAULT_THRESHOLD_FRACTION;
  opt.quantizeBits = DEFAULT_QUANTIZE_BITS;
  opt.quantizeAlgorithm = DEFAULT_QUANTIZE_ALGORITHM;

  for (nextArg = 1; nextArg < argc; nextArg++) {
    const char *arg = argv[nextArg];
    if (arg[0] != '-') break;

    if (!strcmp(arg, "-d")) {
      opt.doCompress = false;
    }

    else if (!strcmp(arg, "-steps")) {
      if (++nextArg >= argc) printHelp();
      arg = argv[nextArg];
      if (1 != sscanf(arg, "%d", &opt.waveletSteps) ||
          opt.waveletSteps < 0) {
        fprintf(stderr, "Invalid # of wavelet transform steps: \"%s\"\n", arg);
        return false;
      }
    }

    else if (!strcmp(arg, "-thresh")) {
      if (++nextArg >= argc) printHelp();
      arg = argv[nextArg];
      if (1 != sscanf(arg, "%f", &opt.thresholdFraction) ||
          opt.thresholdFraction < 0 ||
          opt.thresholdFraction > 1) {
        fprintf(stderr, "Invalid threshold proportion: \"%s\"\n", arg);
        return false;
      }
    }

    else if (!strcmp(arg, "-qbits")) {
      if (++nextArg >= argc) printHelp();
      arg = argv[nextArg];
      if (1 != sscanf(arg, "%d", &opt.quantizeBits) ||
          opt.quantizeBits < 1 ||
          opt.quantizeBits > 32) {
        fprintf(stderr, "Invalid # quantize bits: \"%s\"\n", arg);
        return false;
      }
    }

    else if (!strcmp(arg, "-qalg")) {
      if (++nextArg >= argc) printHelp();
      opt.quantizeAlgorithm = quantAlgName2Id(argv[nextArg]);
      if (opt.quantizeAlgorithm < 0) {
        fprintf(stderr, "Invalid quantize algorithm: \"%s\"\n", argv[nextArg]);
        return false;
      }
    }

    else {
      fprintf(stderr, "Unrecognized option: \"%s\"\n", arg);
      return false;
    }
  }

  return true;
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

  vector<float> quantBinBoundaries;
  vector<float> quantBinValues;

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
    quant_count_cpu_sorted(size*size, sortedAbsData, opt.quantizeBits,
			   threshold, quantBinBoundaries, quantBinValues);
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
  if (!writeQuantDataParamStrings(outputFile, fileData)) return false;
    
  elapsed = NixTimer::time() - startTime;
  printf("Write data file: %.2f ms\n", elapsed*1000);

  elapsed = NixTimer::time() - firstStartTime;
  printf("Total: %.2f ms\n", elapsed*1000);
  

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
  if (!readQuantDataParamStrings(inputFile, f)) return false;

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
  if (opt.quantizeAlgorithm != QUANT_ALG_UNIFORM) {
    printf("%s de-quantization not implemented yet.\n", 
	   quantAlgId2Name(opt.quantizeAlgorithm));
    return false;
  }
  dquant_unif_cpu(size, data, f.quantizeBits, f.threshold, f.quantMaxVal);
  printf("Dequantize: %.2f ms\n", (NixTimer::time() - startTime) * 1000);

  // perform inverse wavelet transform
  elapsed = haar_2d(size, data, true, f.waveletSteps);
  printf("Wavelet inverse transform: %.2f ms\n", elapsed);

  // write the reconstructed data
  startTime = NixTimer::time();
  if (!writeDataFile(outputFile, data, size, size, true)) return false;
  printf("Write file: %.2f ms\n", (NixTimer::time() - startTime) * 1000);

  printf("Total: %.2f ms\n", (NixTimer::time() - firstStartTime) * 1000);

  return true;
}


/**
   Format: See the source.
   width (int, 4 bytes)
   height (int, 4 bytes)
   wavelet steps (int, 4 bytes)
   quantize bits (int, 4 bytes)
   threshold value (float, 4 bytes)
   quantization algorithm enum value (int, 4 bytes) 
   maximum value before quantization (float, 4 bytes)
   all the quantized data (width*height*(# of bits) bits, rounded up to a word)
*/
bool writeQuantDataSimple(const char *filename, FileData &f) {

  FILE *outf = fopen(filename, "wb");
  if (!outf) {
    printf("Error writing to \"%s\"\n", filename);
    return false;
  }

  assert(sizeof f.quantizeAlgorithm == sizeof(int));

  // write parameters to the header
  fwrite(&f.width, sizeof(int), 1, outf);
  fwrite(&f.height, sizeof(int), 1, outf);
  fwrite(&f.waveletSteps, sizeof(int), 1, outf);
  fwrite(&f.quantizeBits, sizeof(int), 1, outf);
  fwrite(&f.threshold, sizeof(float), 1, outf);
  fwrite(&f.quantizeAlgorithm, sizeof(int), 1, outf);
  // XXX we will need something different for Lloyd's algorithm codebook
  fwrite(&f.quantMaxVal, sizeof(float), 1, outf);
    
  bool success = writeQuantData(filename, outf, f);

  fclose(outf);

  return success;
}

/**
   Read the data written by writeQuantData.
*/
bool readQuantDataSimple(const char *filename, FileData &f) {

  FILE *inf = fopen(filename, "rb");
  if (!inf) {
    printf("Cannot read \"%s\"\n", filename);
    return false;
  }

  assert(sizeof f.quantizeAlgorithm == sizeof(int));

  fread(&f.width, sizeof(int), 1, inf);
  fread(&f.height, sizeof(int), 1, inf);
  fread(&f.waveletSteps, sizeof(int), 1, inf);
  fread(&f.quantizeBits, sizeof(int), 1, inf);
  fread(&f.threshold, sizeof(float), 1, inf);
  fread(&f.quantizeAlgorithm, sizeof(int), 1, inf);
  fread(&f.quantMaxVal, sizeof(float), 1, inf);

  bool success = readQuantData(filename, inf, f);

  fclose(inf);

  return success;
}


bool writeQuantDataParamStrings(const char *filename, FileData &f) {

  FILE *outf = fopen(filename, "wb");
  if (!outf) {
    printf("Error writing to \"%s\"\n", filename);
    return false;
  }

  // write parameters to the header
  ParamString p;
  p.setInt("w", f.width);
  p.setInt("h", f.height);
  p.setInt("ws", f.waveletSteps);
  p.setInt("qb", f.quantizeBits);
  p.setFloat("th", f.threshold);
  p.set("qa", quantAlgId2Name(f.quantizeAlgorithm));
  if (f.quantizeAlgorithm == QUANT_ALG_UNIFORM ||
      f.quantizeAlgorithm == QUANT_ALG_LOG) {
    p.setFloat("max", f.quantMaxVal);
  } else {
    p.setFloatList("qbound", f.quantBinBoundaries);
    p.setFloatList("cb", f.quantBinValues);
  }

  p.writeParameters(outf);
    
  bool success = writeQuantData(filename, outf, f);

  fclose(outf);

  return success;
}


bool readQuantDataParamStrings(const char *filename, FileData &f) {

  FILE *inf = fopen(filename, "rb");
  if (!inf) {
    printf("Cannot read \"%s\"\n", filename);
    return false;
  }

  ParamString p;
  p.readParameters(inf);
  if (!p.getInt("w", f.width)) printf("width not found\n");
  if (!p.getInt("h", f.height)) printf("height not found\n");
  if (!p.getInt("ws", f.waveletSteps)) printf("wavelet steps not found\n");
  if (!p.getInt("qb", f.quantizeBits)) printf("quant bits not found\n");
  if (!p.getFloat("th", f.threshold)) printf("threshold not found\n");
  string qa;
  if (!p.get("qa", qa)) printf("quant alg not found\n");
  f.quantizeAlgorithm = quantAlgName2Id(qa.c_str());
  if (f.quantizeAlgorithm == QUANT_ALG_UNIFORM ||
      f.quantizeAlgorithm == QUANT_ALG_LOG) {
    if (!p.getFloat("max", f.quantMaxVal)) printf("max value not found\n");
  }

  bool success = readQuantData(filename, inf, f);

  fclose(inf);

  return success;
}


bool readQuantData(const char *filename, FILE *inf, FileData &f) {
  int count = f.width * f.height;
  f.data = new float[count];

  BitStreamReader bits(inf);
  float *writePos = f.data;
  float *endPos = f.data + count;

  while (writePos < endPos) {
    if (bits.isEmpty()) {
      printf("Ran out of data reading %s, expected %d entries, got %d\n",
             filename, f.width * f.height, (int)(writePos - f.data));
      return false;
    }

    // read sign bit, value, length
    int sign = bits.read(1);
    // map 1,0 -> 2,0 -> -2,0 -> -1,1
    sign = 1 - (sign * 2);

    int quantized = bits.read(f.quantizeBits);
    int value = sign * quantized;

    int length = bits.read(8);
    
    for (int i=0; i < length; i++)
      *writePos++ = (float)value;
  }
  return true;
}


bool writeQuantData(const char *filename, FILE *outf, FileData &f) {

  // write the data
  BitStreamWriter bits(outf);
  int count = f.width*f.height;


  // this object takes (length,value) run-length pairs and writes
  // them in binary to the given bit stream
  WriteRLEPairsToBitStream rleToBits(&bits, f.quantizeBits);

  // this object takes data values as input, and passes (length,value)
  // run-length pairs to the rleToBits object
  EncodeRunLength<WriteRLEPairsToBitStream> rleEncoder(&rleToBits);

  // XXX see if it's faster to use pointer to traverse f.data[]

  for (int i=0; i < count; i++) {
    int x = (int) f.data[i];
    rleEncoder.data(x);
  }
  rleEncoder.end();

  printf("%d input values, %d RLE output pairs\n", count,
  rleEncoder.getOutputCount());
  return true;
}


QuantizeAlgorithm quantAlgName2Id(const char *name) {
  if (!strcmp(name, "uniform")) return QUANT_ALG_UNIFORM;
  if (!strcmp(name, "log")) return QUANT_ALG_LOG;
  if (!strcmp(name, "lloyd")) return QUANT_ALG_LLOYD;
  if (!strcmp(name, "count")) return QUANT_ALG_COUNT;
  return QUANT_ALG_UNKNOWN;
}
  

const char *quantAlgId2Name(QuantizeAlgorithm id) {
  switch (id) {
  case QUANT_ALG_UNIFORM: return "uniform";
  case QUANT_ALG_LOG: return "log";
  case QUANT_ALG_COUNT: return "count";
  case QUANT_ALG_LLOYD: return "lloyd";
  default: return NULL;
  }
}
