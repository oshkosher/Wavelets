#ifndef __TEST_COMPRESS_COMMON_H__
#define __TEST_COMPRESS_COMMON_H__

#include <vector>
#include <string>
#include "wavelet.h"
#include "dwt_cpu.h"
#include "data_io.h"
#include "bit_stream.h"
#include "wavelet.h"
#include "wavelet_compress.pb.h"
#include "cubelet_file.h"

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



void printHelp();
bool parseOptions(int argc, char **argv, Options &opt, int &nextArg);

// Write 'cube' to a cubelet stream.
// If sizeBytes is not NULL, store the size of the output data in it.
bool writeQuantData(CubeletStreamWriter &cubeletStream,
                    CubeInt *cube, Options &opt,
                    int *sizeBytes = NULL);

// Read cubelet from a file
bool readQuantData(CubeletStreamReader &cubeletStream, CubeInt *data);


#endif // __TEST_COMPRESS_COMMON_H__
