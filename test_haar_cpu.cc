#include <cstdio>
#include "dwt_cpu.h"
#include "data_io.h"


int main3() {
  float data[] = {9, 7, 3, 5, 6, 10, 2, 6};
  int len = 8;

  print_matrix(len, 1, data);
  printf("\n");

  haar_not_lifting(len, data, false, 1);
  print_matrix(len, 1, data);

  printf("\n");
  haar_not_lifting(len, data, true, 1);
  print_matrix(len, 1, data);

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
  

int main(int argc, char **argv) {
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

  if (width != height) {
    printf("Only square matrices are supported.\n");
    return 1;
  }

  int size = height;

  if (stepCount < 0) 
    haar_not_lifting_2d(size, data, true, -stepCount);
  else
    haar_not_lifting_2d(size, data, false, stepCount);
  
  // printMatrix(width, height, data);
  if (writeDataFile(outputFile, data, size, size, 0))
    printf("%s written\n", outputFile);
  
  return 0;
}
