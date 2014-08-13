#include <cstdio>
#include <cassert>
#include "dwt_cpu.h"
#include "data_io.h"


int testPad() {
  /*
  // float data[] = {9, 7, 3, 5, 6, 10, 2, 6};
  float *output = NULL, data[] = {17, 15, 13, 20, 35};
  int len = 5, outputLen;

  print_matrix(len, 1, data);
  printf("\n");

  outputLen = dwt_padded_length(len, 4, true);
  output = dwt_pad(len, data, outputLen, output, REFLECT);
  print_matrix(outputLen, 1, output);

  printf("\n");
  haar_not_lifting(len, output, true, 1);
  print_matrix(len, 1, output);
  */

  float data[] = {
    17,   14,    1,   16,   20,    6,    2,   11,
    14,    9,    4,   10,    6,   16,   17,   13,
     3,    7,   17,   11,   10,    5,    2,   16,
     6,   12,   15,   13,   19,    7,   13,   20,
    13,   19,    2,   13,    6,    8,   10,    8,
    10,    8,   13,   15,    2,   15,   12,   11,
    20,   13,   17,   20,   20,    2,   14,    3,
     8,   12,    2,   19,   15,   14,    8,   16 
  };

  dwt_pad_2d(5, 5, 8, data, 8, 8, 8, data, REFLECT);
  print_matrix(8, 8, data);

  return 0;
}



int main2() {

  // float data[] = {9, 7, 3, 5, 6, 10, 2, 6};
  float data[] = {
    17,   14,    1,   16,   20,    6,    2,   11,
    14,    9,    4,   10,    6,   16,   17,   13,
     3,    7,   17,   11,   10,    5,    2,   16,
     6,   12,   15,   13,   19,    7,   13,   20,
    13,   19,    2,   13,    6,    8,   10,    8,
    10,    8,   13,   15,    2,   15,   12,   11,
    20,   13,   17,   20,   20,    2,   14,    3,
     8,   12,    2,   19,   15,   14,    8,   16 
  };
  int len = 8;


  /*
  haar_not_lifting(len, data);
  printf("\n");
  haar_inv_not_lifting(len, data);
  */

  haar_not_lifting_2d(len, data, false, 1);
  print_matrix(len, len, data);

  printf("\n");

  haar_not_lifting_2d(len, data, true, 1);
  print_matrix(len, len, data);

  return 0;
}
  

int main_full(int argc, char **argv) {
  if (argc != 4) {
    fprintf(stderr, "\n  test_haar_cpu <steps> <input data> <output data>\n"
            "  Negative steps does inverse transform\n\n");
    return 1;
  }

  float *data;
  int width, height, stepCount;
  if (1 != sscanf(argv[1], "%d", &stepCount)) {
    printf("Invalid step count: \"%s\"\n", argv[1]);
    return 1;
  }
  const char *inputFile = argv[2], *outputFile = argv[3];
  printf("Reading %s...", inputFile);
  fflush(stdout);
  if (!readDataFile(inputFile, &data, &width, &height)) return 1;
  printf("%d x %d\n", width, height);
  fflush(stdout);

  // pad the size to a power of two
  int size = (width > height) ? width : height;
  int padSize = dwt_padded_length(size, 0, true);
  if (width != size || height != size) {
    printf("Padding to %dx%d\n", padSize, padSize);
    fflush(stdout);
    float *padData = dwt_pad_2d(height, width, width, data,
				padSize, padSize, padSize, NULL, REPEAT);
    delete[] data;
    data = padData;
    size = width = height = padSize;
  }

  if (width != height) {
    printf("Only square matrices are supported.\n");
    return 1;
  }

  if (stepCount < 0) 
    haar_not_lifting_2d(size, data, true, -stepCount);
  else
    haar_not_lifting_2d(size, data, false, stepCount);
  
  // printMatrix(width, height, data);
  printf("Writing...\n");
  fflush(stdout);
  if (writeDataFile(outputFile, data, size, size, true))
    printf("%s written\n", outputFile);

  delete[] data;
  
  return 0;
}


void testMisc() {
  // test countLeadingZeros
  assert(countLeadingZeros(0) == 32);
  assert(countLeadingZeros(1) == 31);
  assert(countLeadingZeros(2) == 30);
  assert(countLeadingZeros(3) == 30);
  assert(countLeadingZeros(0x7fffffff) == 1);
  assert(countLeadingZeros(0xffffffff) == 0);

  // test ceilLog2()
  assert(ceilLog2(1) == 0);
  assert(ceilLog2(2) == 1);
  assert(ceilLog2(3) == 2);
  assert(ceilLog2(4) == 2);
  assert(ceilLog2(5) == 3);
  assert(ceilLog2(6) == 3);
  assert(ceilLog2(7) == 3);
  assert(ceilLog2(8) == 3);
  assert(ceilLog2(9) == 4);
  assert(ceilLog2(511) == 9);
  assert(ceilLog2(512) == 9);
  assert(ceilLog2(513) == 10);

  // test dwt_padded_length()
  assert(dwt_padded_length(4, 1, false) == 4);
  assert(dwt_padded_length(5, 1, false) == 6);
  assert(dwt_padded_length(5, 2, false) == 8);
  assert(dwt_padded_length(6, 1, false) == 6);
  assert(dwt_padded_length(6, 2, false) == 8);
  assert(dwt_padded_length(3, 0, true) == 4);
  assert(dwt_padded_length(4, 0, true) == 4);
  assert(dwt_padded_length(5, 0, true) == 8);

  printf("OK.\n");
}

int main(int argc, char **argv) {
  // testMisc();
  // testPad();

  main_full(argc, argv);

  return 0;
}
