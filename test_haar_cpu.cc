#include <cstdio>
#include "dwt_cpu.h"
#include "data_io.h"


int main2() {

  float data[] = {9, 7, 3, 5, 6, 10, 2, 6};
  int len = sizeof data / sizeof(float);

  print_matrix(len, 1, data);

  haar_not_lifting(len, data);

  printf("\n");

  haar_inv_not_lifting(len, data);

  return 0;
}
  

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "\n  test_haar_cpu <input data> <output data>\n\n");
    return 1;
  }

  float *data;
  int width, height;
  printf("Reading %s...", argv[1]);
  fflush(stdout);
  if (!readDataFile(argv[1], &data, &width, &height)) return 1;
  printf("%d x %d\n", width, height);
  fflush(stdout);

  haar_not_lifting_2d(height, width, data, true);
  
  // printMatrix(width, height, data);
  if (writeDataFile(argv[2], data, width, height, 1))
    printf("%s written\n", argv[2]);
  
  
  return 0;
}
