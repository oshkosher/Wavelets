#ifndef __TEST_COMPRESS_COMMON_H__
#define __TEST_COMPRESS_COMMON_H__

#include <cmath>
#include <vector>
#include <string>
#include "wavelet.h"
#include "dwt_cpu.h"
#include "data_io.h"
#include "bit_stream.h"
#include "wavelet.h"
#include "wavelet_compress.pb.h"
#include "cubelet_file.h"
#include "quant.h"
#include "optimize.h"

#define FILE_ID_STRING "SCU wavelet 1.0\n"

// global variable that disables status output
extern bool QUIET;

/** This holds the the parameters as set by the user on the command line,
    and the ones saved in the data file. */
struct Options {

  // if false, then decompress
  bool doCompress;

  // suppress output
  bool quiet;

  // greatly expand output
  bool verbose;

  // After computing the bit encoding of each value, print the encodings
  bool printHuffmanEncoding;

  // After quantizing, dequant, reverse the wavelet transform, and compare
  // the original image with the compressed & decompressed result
  bool doComputeError;

  // After loading the data and doing the wavelet transform, call
  // the parameter optimzation routine to set the threshold and bin count.
  bool doOptimize;

  // -bq option : if not "", save a copy data in this file before quantizing
  std::string saveBeforeQuantizingFilename;

  bool runQuantizationExperiments;

  // all the parameters for the wavelet compression
  WaveletCompressionParam param;

  Options() {
    init();
  }

  void init() {
    doCompress = true;
    quiet = false;
    verbose = false;
    printHuffmanEncoding = false;
    doComputeError = false;
    doOptimize = false;
    saveBeforeQuantizingFilename = "";
    runQuantizationExperiments = false;

    param.init();
  }
};


/** Feed this pairs of data (version1, version2), and it will
    produce summaries of the differences.
*/
class ErrorAccumulator {

  // declare both these as double to reduce the accumulated error
  double sumDiff;         // sum of the absolute differences
  double sumDiffSquared;  // sum of the squares of the absolute differences

  float maxValue;
  int count;

 public:
 ErrorAccumulator() : sumDiff(0), sumDiffSquared(0), maxValue(0), count(0) {}

  // when computing peak signal-to-noise ratio, we need to know the
  // maximum possible value.
  void setMaxPossible(float m) {
    if (m <= 0) {
      fprintf(stderr, "ERROR: invalid maximum value (%f) for "
              "peak signal-to-noise ratio calculation.\n", m);
    } else {
      maxValue = m;
    }
  }

  void setSumDiff(double s) {sumDiff = s;}
  void setSumDiffSquared(double s) {sumDiffSquared = s;}
  void setCount(int c) {count = c;}

  // add one pair of values
  void add(float a, float b) {
    float diff = fabsf(a-b);
    sumDiff += diff;
    sumDiffSquared += diff*diff;
    count++;
  }

  int getCount() {return count;}

  float getAverageError() {
    if (count == 0) return 0;
    return (float)(sumDiff / count);
  }

  float getL1Error() {return (float)sumDiff;}
  
  float getL2Error() {return (float)sqrt(sumDiffSquared);}
  
  float getMeanSquaredError() {
    if (count == 0) return 0;
    return (float)(sumDiffSquared / count);
  }

  float getPeakSignalToNoiseRatio() {
    if (maxValue == 0) return 0;
    float mse = getMeanSquaredError();
    if (mse == 0) return 0;
    return 10 * log10(maxValue * maxValue / mse);
  }
};


void printHelp();
bool parseOptions(int argc, char **argv, Options &opt, int &nextArg);

// returns true iff 'suffix' is a suffix of 's'
bool endsWith(const char *s, const char *suffix);

// Read one cubelet from the input file into 'inputData' (the original data)
// and a copy in 'data' (padded and translated to floats).
// Both .data (data_io.h) and .cube (cubelet_file.h) files are supported.
// Based on the compression parameters, the data may be padded also.
bool readData(const char *filename,
              Cube *inputData,  // the original data, may not be floats
              CubeFloat *data,  // input data transformed into floats
              CubeletStreamReader *cubeletStream,
              WaveletCompressionParam *param);

// If the user specified too many steps or a negative number of steps,
// set them to the maximum number given the data size.
void setWaveletSteps(scu_wavelet::int3 &steps, const scu_wavelet::int3 &size);

// Find the size to which the data needs to be padded so each transform
// step has an even number of elements. Change the given size in-place.
void padDataSize(scu_wavelet::int3 &size, scu_wavelet::int3 steps);

// pad the data, repeating the last value in each axis as needed
void padData(CubeFloat &data, scu_wavelet::int3 originalSize);

// make the cubelet data into floats, if it isn't already
void translateCubeDataToFloat(Cube *src, CubeFloat *dest);


// Peform the wavelet transform
bool waveletTransform(CubeFloat &data, const WaveletCompressionParam &param,
                      bool isInverse, bool verbose);

bool dequantize(const CubeInt &quantizedData, CubeFloat &data,
                const WaveletCompressionParam &param);

void computeErrorRates(const CubeInt *quantizedData,
                       const WaveletCompressionParam &param,
                       const Cube *inputData,
                       float *meanSquaredError,
                       float *peakSNR);


// this will modify restoredData in place
void computeErrorRatesAfterDequant
(CubeFloat *restoredData,
 const WaveletCompressionParam &param,
 const Cube *inputData,
 ErrorAccumulator *errAccum);

// turn the data back into the original datatype, if it wasn't floats
void translateCubeDataToOriginal(CubeFloat *src, Cube *dest,
                                 bool verbose = false);

// Read cubelet from a file
bool readQuantData(CubeletStreamReader &cubeletStream, CubeInt *data);

// Write 'cube' to a cubelet stream.
// If sizeBytes is not NULL, store the size of the output data in it.
bool writeQuantData(CubeletStreamWriter &cubeletStream,
                    CubeInt *cube, Options &opt,
                    int *sizeBytes = NULL);

Quantizer *createQuantizer(const WaveletCompressionParam &param);

#endif // __TEST_COMPRESS_COMMON_H__
