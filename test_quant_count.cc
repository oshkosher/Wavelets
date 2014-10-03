#include <cstdio>
#include <vector>
#include "quant_count.h"
#include "data_io.h"


int main2() {
  float data[100];
  int len = 100;
  for (int i=0; i < 100; i++) data[i] = i / 100.0f;

  std::vector<float> boundaries, codebook;

  quant_count_cpu(len, data, 2, .4, boundaries, codebook);
  for (int i=0; i < 3; i++)
    printf("  %f", boundaries[i]);
  printf("\n");
  for (int i=0; i < 4; i++)
    printf("  %f", codebook[i]);
  printf("\n");

  return 0;
}


int main() {
  float *data;
  int width, height;
  
  if (!readDataFile("lenna.data", &data, &width, &height)) {
    printf("Failed to read lenna.data\n");
    return 1;
  }

  std::vector<float> boundaries, codebook;
  int bits = 3;

  quant_count_cpu(width * height, data, bits, .1, boundaries, codebook);
  int binCount = 1 << bits;
  for (int i=0; i < binCount-1; i++)
    printf("  %f", boundaries[i]);
  printf("\n");
  for (int i=0; i < binCount; i++)
    printf("  %f", codebook[i]);
  printf("\n");

  return 0;
}
