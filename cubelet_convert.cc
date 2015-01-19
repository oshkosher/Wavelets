#include <string>
#include <cstdio>
#include <cstdlib>
#include <cinttypes>
#include "cubelet_file.h"

using namespace std;
using namespace scu_wavelet;

/*
  Convert between our cubelet format and other formats.


  List all the cubelets in a cubelet file.
    cubelet_convert foo.cube

  Convert a file of unsigned bytes (values 0..255) into a cubelet file:
    cubelet_convert -bytes data.img 512 128 1024 foo.cube

  Slice up the input:
    cubelet_convert -bytes 512,128,1024 -cubesize 128,128,128 data.img foo.cube

  Slice with overlapping cubelets
    cubelet_convert -bytes 512,128,1024 -cubesize 128,128,128 -overlap 4,4,4 \
      data.img foo.cube
    
  Convert the first cubelet to raw data
    cubelet_convert foo.cube -bytes data.img

  Convert a different cubelet to raw data
    cubelet_convert -cubeid 1,2,3 foo.cube -bytes data.img


  The values 512, 128, and 1024 are the width, height, and depth
  of the given data.
*/

struct Int3 {
  int x, y, z;

  // given a string in the form "123,456,789" parse it into up to
  // three integers, store them into x, y, and z, and return the number
  // of values successfully parsed. Set unparsed values to 1.
  int parse(const char *str) {
    x = y = z = 1;
    return sscanf(str, "%d,%d,%d", &x, &y, &z);
  }

  void set(int x_, int y_, int z_) {
    x = x_;
    y = y_;
    z = z_;
  }

  bool isValid() {
    return x > 0 && y > 0 && z > 0;
  }
};

struct Options {
  std::string inputFile;
  int3 inputShape;

  std::string outputFile;
};
  

void printHelp();
bool parseOptions(Options &opt, int argc, const char **argv);
bool listCubeletFile(const std::string &filename);
bool isCubeletFilename(const std::string &filename);
bool convertFromRawBytes(Options &opt);


int main(int argc, const char **argv) {
  Options opt;

  if (!parseOptions(opt, argc, argv)) return 1;

  if (opt.outputFile.empty()) {
    listCubeletFile(opt.inputFile);
  }

  else {
    if (isCubeletFilename(opt.outputFile) &&
        !isCubeletFilename(opt.inputFile)) {
      convertFromRawBytes(opt);
    }
  }
  
  return 0;
}


void printHelp() {
  printf("\n  cubelet_convert [input opt] <inputfile> [output opt] <outputfile>\n"
         "    to do file conversion\n"
         "  cubelet_convert <.cube_input_file>\n"
         "    to list cubelets in a flie\n"
         "\n  Input options:\n"
         "  -bytes <w,h,d>\n"
         "\n");
  exit(1);
}

bool parseOptions(Options &opt, int argc, const char **argv) {
  opt.inputFile = "";
  opt.outputFile = "";
  opt.inputShape = int3(-1,-1,-1);

  bool inputFileRead = false;
  
  for (int argno = 1; argno < argc; argno++) {

    const char *arg = argv[argno];

    if (!strcmp(arg, "-bytes")) {
      arg = argv[++argno];
      if (!arg) printHelp();

      if (inputFileRead) {
        fprintf(stderr, "Error: only the shape of the input file can be "
                "specified.\n");
        return false;
      }

      int count = opt.inputShape.parse(arg);
      if (count < 1) {
        fprintf(stderr, "Couldn't parse input file shape \"%s\". "
                "Expecting X,Y,Z.\n", arg);
        return false;
      }
    }
    
    else {  // just a filename

      if (!inputFileRead) {
        opt.inputFile = arg;
        inputFileRead = true;
      } else {
        opt.outputFile = arg;
      }
    }

  }

  if (!inputFileRead) printHelp();
  
  return true;
}


// return true iff the given filename has a ".cube" suffix.
bool isCubeletFilename(const string &filename) {
  return filename.length() >= 5 &&
    filename.substr(filename.length() - 5, 5) == ".cube";
}


bool listCubeletFile(const string &filename) {

  CubeletStreamReader in;

  if (!in.open(filename.c_str())) return false;

  Cube cube;
  
  while (true) {
    if (!in.next(&cube)) break;

    const char *compressed = cube.isWaveletCompressed ? " compressed" : "";

    printf("Cubelet %d,%d,%d: %dx%dx%d of %s, %d%s bytes at %" PRIu64 "\n",
           cube.parentOffset.x, cube.parentOffset.y, cube.parentOffset.z,
           cube.width(), cube.height(), cube.depth(),
           waveletDataTypeName(cube.datatype),
           cube.getSizeInBytes(), compressed, cube.dataFileOffset);
  }

  return true;
}


bool convertFromRawBytes(Options &opt) {

  if (!(opt.inputShape > int3(0,0,0))) {
    fprintf(stderr, "Input shape is unspecified or invalid.\n");
    return false;
  }
  
  FILE *inf = fopen(opt.inputFile.c_str(), "rb");
  if (!inf) {
    fprintf(stderr, "Failed to read %s\n", opt.inputFile.c_str());
    return false;
  }

  CubeByte cube;
  cube.size = int3(opt.inputShape.x, opt.inputShape.y, opt.inputShape.z);
  cube.allocate();
  uint64_t dataSize = (uint64_t)opt.inputShape.x * opt.inputShape.y
    * opt.inputShape.z;
  uint64_t bytesRead;
  bytesRead = fread(cube.data(), 1, dataSize, inf);
  if (bytesRead != dataSize) {
    fprintf(stderr, "Data size mismatch: expected %" PRIu64 " bytes, "
            "got %" PRIu64 "\n", dataSize, bytesRead);
  }

  CubeletStreamWriter out;
  if (!out.open(opt.outputFile.c_str())) return false;
  if (!out.addCubelet(&cube)) return false;
  out.close();

  printf("Wrote %s\n", opt.outputFile.c_str());

  return true;
}

