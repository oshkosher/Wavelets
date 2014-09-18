/*
  Read a data file and output all the data values in a single column,
  to make further processing easy.
*/

#include <cstdio>
#include <cstdlib>
#include "data_io.h"

void printHelp() {
  printf("\n  list_data <inputfile>\n\n");
  exit(1);
}

int main(int argc, char **argv) {
  if (argc != 2) printHelp();

  float *data;
  int width, height;
  if (!readDataFile(argv[1], &data, &width, &height)) {
    printf("Failed to read %s.\n", argv[1]);
    return 1;
  }

  int count = width * height;
  for (int i=0; i < count; i++) {
    printf("%g\n", data[i]);
  }

  return 0;
}
