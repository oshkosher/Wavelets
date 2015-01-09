#include <cstdio>
#include "data_io.h"

/*
  Normalize recursive quadrants of an image, to make a nicer visualization
  of an image process by a discrete wavelet transform.

  Take a 1-step transform, where the results are four quadrants:
    A  B
    C  D
  An approximmation of the image will be in quadrant A.
  B, C, and D will contain differential values.
  The pixels in A will be similar to the image. If the image consisted of
  values from 0 to 1, these will be in the range 0..2. However, B, C, and D
  may be in the range -1..1.

  This tool processes each of the quadrants, doing a linear remapping
  of the values to the range 0..1.

  In a multi-step transform, the quadrants will be processed recursively:

    A  B0  C0 C0
    B1 B2  C0 C0

    C1 C1  C2 C2
    C1 C1  C2 C2

 */


inline float getPixel(int x, int y, float data[], int width) {
  return data[y*width + x];
}

inline void setPixel(float value, int x, int y, float data[], int width) {
  data[y*width + x] = value;
}

void normalizeRegion(float data[], int totalWidth,
                     int x, int width,
                     int y, int height) {
  
  float min = getPixel(x, y, data, totalWidth);
  float max = min;
  
  for (int j=0; j < height; j++) {
    for (int i=0; i < width; i++) {
      float p = getPixel(x + i, y + j, data, totalWidth);
      if (p < min) min = p;
      if (p > max) max = p;
    }
  }

  printf("normalize at (%d,%d) size (%d,%d)  %g .. %g\n",
         x, y, width, height, min, max);

  float offset = -min;
  float scale = (max == min) ? 0 : 1 / (max-min);

  for (int j=0; j < height; j++) {
    for (int i=0; i < width; i++) {
      float p = getPixel(x + i, y + j, data, totalWidth);
      setPixel((p + offset) * scale - .5f, x+i, y+j, data, totalWidth);
    }
  }
}  


/*
  Initial call: steps = 0, step = 0
  One call, the whole region

  Initial call: steps = 1, step = 0
    a b
    e f
    recursive call 1, 1: process a
    then the three regions: b, e, f

  Initial call: steps = 2, step = 0
    a b c d
    e f g h
    i j k l
    m n o p

    recursive call 2, 1: process (abef)
    then the three regions: (cdgh), (ijmn), (klop)
 */
void recursiveNormalize(float data[], int width, int height,
                         int steps, int step) {

  if (steps == step) {
    normalizeRegion(data, width, 0, width>>step, 0, height>>step);
    return;
  }

  recursiveNormalize(data, width, height, steps, step+1);
  
  int myWidth = width >> (step+1);
  int myHeight = height >> (step+1);

  // right
  normalizeRegion(data, width, myWidth, myWidth, 0, myHeight);
  // down
  normalizeRegion(data, width, 0, myWidth, myHeight, myHeight);
  // right/down
  normalizeRegion(data, width, myWidth, myWidth, myHeight, myHeight);
}


  



int main(int argc, char **argv) {
  const char *inputFile, *outputFile;
  int steps, width, height;
  float *data;

  if (argc != 4) {
    printf("\n  normalize <steps> <inputdata> <outputdata>\n\n");
    return 1;
  }

  if (1 != sscanf(argv[1], "%d", &steps)) {
    printf("Invalid values for steps: \"%s\"\n", argv[1]);
    return 1;
  }

  inputFile = argv[2];
  outputFile = argv[3];

  if (!readDataFile(inputFile, &data, &width, &height)) return 1;

  recursiveNormalize(data, width, height, steps, 0);
  
  writeDataFile(outputFile, data, width, height, true);

  return 0;
}
