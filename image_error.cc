#include <cstdio>
#include <cstdlib>
#include "data_io.h"

/*
  Compute the mean squared error between two images.
  They must be the same size.
*/

void printHelp();
#define MIN(a,b) ((a) <= (b) ? (a) : (b))


int main(int argc, char **argv) {
  
  if (argc != 3) printHelp();
  
  double *data1, *data2;
  int width1, width2, height1, height2;
  const char *file1 = argv[1], *file2 = argv[2];

  if (!readDataFile(argv[1], &data1, &width1, &height1)) return 1;
  if (!readDataFile(argv[2], &data2, &width2, &height2)) return 1;

  int width = MIN(width1, width2), height = MIN(height1, height2);
  if (width1 != width2 || height1 != height2) {
    printf("Image size mismatch. %s %dx%d, %s %dx%d.\n"
           "Comparing just the upper %dx%d pixels.\n",
           file1, height1, width1,
           file2, height2, width2,
           width, height);
  }

  double sumSquaredErr = 0;
  
  for (int y=0; y < height; y++) {
    for (int x=0; x < width; x++) {
      double err = data1[y*width1 + x] - data2[y*width2 + x];
      sumSquaredErr += err*err;
    }
  }

  delete[] data1;
  delete[] data2;

  double meanSquaredErr = sumSquaredErr / ((double)width * height);

  printf("Mean squared error: %.5g\n", meanSquaredErr);

  return 0;
}


void printHelp() {
  fprintf(stderr, "\n  image_error <file1> <file2\n"
          "  Prints the mean squared error between the two data files.\n\n");
  exit(1);
}
