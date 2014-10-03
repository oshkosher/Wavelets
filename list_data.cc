/*
  Read a data file and output all the data values in a single column,
  to make further processing easy.
*/

#include <cstdio>
#include <cstdlib>
#include "data_io.h"

void printHelp() {
  printf("\n  list_data <inputfile> [<xOffset> <yOffset> <width> <height>]\n\n");
  exit(1);
}

int main(int argc, char **argv) {
  if (argc != 2 && argc != 6) printHelp();

  float *data;
  int width, height;
  if (!readDataFile(argv[1], &data, &width, &height)) {
    fprintf(stderr, "Failed to read %s.\n", argv[1]);
    return 1;
  }

  int xOffset=0, yOffset=0, extractWidth=width, extractHeight=height;
  if (argc == 6) {
    if (1 != sscanf(argv[2], "%d", &xOffset) ||
	1 != sscanf(argv[3], "%d", &yOffset) ||
	1 != sscanf(argv[4], "%d", &extractWidth) ||
	1 != sscanf(argv[5], "%d", &extractHeight)) printHelp();
    if (xOffset < 0 || yOffset < 0 || extractWidth < 0 || extractHeight < 0 ||
	xOffset + extractWidth > width ||
	yOffset + extractHeight > height) {
      fprintf(stderr, "Invalid range for %dx%d image.\n", width, height);
      return 1;
    }
  }

  for (int y=yOffset; y < yOffset + extractHeight; y++) {
    int pos = y * width + xOffset;
    for (int i=0; i < extractWidth; i++) {
      printf("%g\n", data[pos++]);
    }
  }
  /*
  int count = width * height;
  for (int i=0; i < count; i++) {
    printf("%g\n", data[i]);
  }
  */
  return 0;
}
