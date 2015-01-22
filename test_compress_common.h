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

#define FILE_ID_STRING "SCU wavelet 1.0\n"


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
    saveBeforeQuantizingFilename = "";
    runQuantizationExperiments = false;

    param.init();
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
    return sumDiff / count;
  }

  float getL1Error() {return sumDiff;}
  
  float getL2Error() {return sqrt(sumDiffSquared);}
  
  float getMeanSquaredError() {
    if (count == 0) return 0;
    return sumDiffSquared / count;
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

// Write 'cube' to a cubelet stream.
// If sizeBytes is not NULL, store the size of the output data in it.
bool writeQuantData(CubeletStreamWriter &cubeletStream,
                    CubeInt *cube, Options &opt,
                    int *sizeBytes = NULL);

// Read cubelet from a file
bool readQuantData(CubeletStreamReader &cubeletStream, CubeInt *data);

Quantizer *createQuantizer(const WaveletCompressionParam &param);

#endif // __TEST_COMPRESS_COMMON_H__
