#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include "wavelet.h"
#include "cubelet_file.h"

using namespace scu_wavelet;

/*
  This is a one-off program used to conver the "Photron" series of images.

  The input is a series of 16-bit raw data files, and it's turned into
  a cubelet file where each cubelet is one frame of the video, so a
  2-d cubelet with depth==1.

  The input data ranges from 0..65535, but only the values (0 .. 4095)
  * 16 are used, so the data is squished down to values in the range
  0..4095.

  Also, a couple frames were incomplete, so for incomplete frames, fill the
  remaining pixels with 0 so all the frames in the cubelet file are identical.
*/

// input values are (0 .. 4095) * 16
const int MAX_POSSIBLE_VALUE = 4095;
int convertInputValue(int x) {
  return x >> 4;
}

void printHelp() {
  printf("\n  shortImagesToCubelets <width>,<height> <output.cube> [input files...]\n\n");
  exit(1);
}


int main(int argc, char **argv) {
  if (argc < 4) printHelp();

  CubeInt frame;

  if (2 != frame.size.parse(argv[1])) printHelp();
  if (frame.size.x < 1 || frame.size.y < 1) {
    printf("Invalid size: %dx%d\n", frame.size.x, frame.size.y);
    return 1;
  }
  frame.size.z = 1;
  if (!frame.allocate()) return 1;
  int pixelsPerFrame = frame.size.count();
  frame.maxPossibleValue = MAX_POSSIBLE_VALUE;
  
  const char *outputCubeletFilename = argv[2];
  
  CubeletStreamWriter cubeStream;
  if (!cubeStream.open(outputCubeletFilename)) return 1;

  uint16_t *inputFrame = new uint16_t[pixelsPerFrame];
  if (!inputFrame) {
    fprintf(stderr, "Out of memory\n");
    return 1;
  }

  int frameCount = argc - 3;
  for (int argno = 3; argno < argc; argno++) {
    int frameNo = argno - 3;

    // clear the buffer
    memset(inputFrame, 0, 2 * pixelsPerFrame);

    const char *inputFile = argv[argno];
    FILE *inf = fopen(inputFile, "rb");
    if (!inf) {
      printf("Failed to open %s\n", inputFile);
      continue;
    }

    // read the data
    // if the file is too short, just output a warning
    int pixelsRead = fread(inputFrame, 2, pixelsPerFrame, inf);
    if (pixelsRead < pixelsPerFrame) {
      printf("%s too small: %d of %d pixels read\n",
             inputFile, pixelsRead, pixelsPerFrame);
    }
    fclose(inf);

    uint16_t *readPtr = inputFrame;
    int *writePtr = frame.pointer(0,0,0);
    for (int i=0; i < pixelsPerFrame; i++) {
      *writePtr++ = convertInputValue(*readPtr++);
    }

    frame.parentOffset.z = frameNo;
    if (!cubeStream.addCubelet(&frame)) break;
    printf("Frame %d of %d\r", frameNo+1, frameCount);
    fflush(stdout);
  }
  
  if (!cubeStream.close()) return 1;

  printf("%d cubes written to %s\n", frameCount, outputCubeletFilename);

  return 0;
}

  
