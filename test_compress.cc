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
#include "data_io.h"
#include "dwt_cpu.h"

using namespace std;

#define DEFAULT_WAVELET_STEPS 3
#define DEFAULT_THRESHOLD_FRACTION .1
#define DEFAULT_QUANTIZE_BITS 8
#define DEFAULT_QUANTIZE_ALGORITHM "uniform"


struct Options {
  int waveletSteps;
  float thresholdFraction;
  int quantizeBits;
  string quantizeAlgorithm;
};


void printHelp() {
  printf("\n"
         "  test_compress [opts] <input file> <output file>\n"
         "  Options:\n"
         "    -steps <n> : # of wavelet transform steps (default = %d)\n"
         "    -thresh <n> : proportion of data to be eliminated in threshold cutoff\n"
         "                  Must be between [0..1]. (default = %.3f)\n"
         "    -qbits <n> : # of bits into which data is quantized. Must be between \n"
         "                 1 and 32.  (default=%d)\n"
         "    -qalg <alorithm> : quantization algorithm: uniform, log, or lloyd\n"
         "                       (default = %s)\n"
         "\n",
         DEFAULT_WAVELET_STEPS,
         DEFAULT_THRESHOLD_FRACTION,
         DEFAULT_QUANTIZE_BITS,
         DEFAULT_QUANTIZE_ALGORITHM
         );
  exit(1);
}


bool parseOptions(int argc, char **argv, Options &opt, int &nextArg) {
  opt.waveletSteps = DEFAULT_WAVELET_STEPS;
  opt.thresholdFraction = DEFAULT_THRESHOLD_FRACTION;
  opt.quantizeBits = DEFAULT_QUANTIZE_BITS;
  opt.quantizeAlgorithm = DEFAULT_QUANTIZE_ALGORITHM;

  for (nextArg = 1; nextArg < argc; nextArg++) {
    const char *arg = argv[nextArg];
    if (arg[0] != '-') break;

    if (!strcmp(arg, "-steps")) {
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
      opt.quantizeAlgorithm = argv[nextArg];
      if (opt.quantizeAlgorithm != "uniform" &&
          opt.quantizeAlgorithm != "log" &&
          opt.quantizeAlgorithm != "lloyd") {
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

  if (!readDataFile(inputFile, &data, &width, &height)) return 1;

  // pad the data to make it a square power of two, if necessary
  int longSide = width > height ? width : height;
  int size = dwt_padded_length(longSide, 0, true);
  if (width != size || height != size) {
    float *paddedData = dwt_pad_2d(height, width, width, data,
                                   size, size, size,
                                   NULL, REFLECT);
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

  // find the threshold value
  float maxVal, minVal;
  float threshold = thresh_cpu(size, data, opt.thresholdFraction,
                               &maxVal, &minVal);

  // quantize the data
  if (opt.quantizeAlgorithm == "uniform") {
    quant_unif_cpu(size, data, opt.quantizeBits, threshold, maxVal);
  } else if (opt.quantizeAlgorithm == "log") {
    quant_log_cpu(size, data, opt.quantizeBits, threshold, maxVal);
  } else if (opt.quantizeAlgorithm == "lloyd") {
    fprintf(stderr, "Lloyd's algorithm not integrated yet.\n");
    return false;
  } else {
    fprintf(stderr, "Quantization algorithm \"%s\" not found.\n",
            opt.quantizeAlgorithm.c_str());
    return false;
  }

  // write the quantized data to a file
  
}


int main(int argc, char **argv) {

  Options opt;
  int nextArg;

  if (!parseOptions(argc, argv, opt, nextArg)) return 1;

  if (argc - nextArg != 2) printHelp();

  const char *inputFile = argv[nextArg++];
  const char *outputFile = argv[nextArg++];

  if (!processFile(inputFile, outputFile, opt)) return 1;

  return 0;
}
