#include <cstdio>
#include "lloyds.h"


int main() {
  
  float points[] = {1, 2, 3, 10, 11, 12, 20, 21, 22,
                    30, 31, 32, 40, 41, 42, 50, 51, 52};
  int psize = sizeof points / sizeof(float);
  float codebook[] = {10, 20, 30};
  int csize = sizeof codebook / sizeof(float);
  float *partitions = new float[csize-1];
  float dist, reldist;
  unsigned int *quantizedPoints = new unsigned int[psize];
  
  lloyd(points, psize, codebook, csize, partitions, dist, reldist,
        quantizedPoints);

  printf("Codebook: ");
  for (int i=0; i < csize; i++) printf(" %f", codebook[i]);
  putchar('\n');

  printf("Partitions: ");
  for (int i=0; i < csize-1; i++) printf(" %f", partitions[i]);
  putchar('\n');

  printf("Quantized values:\n");
  for (int i=0; i < psize; i++) {
    printf("%f -> %d\n", points[i], quantizedPoints[i]);
  }

  delete[] partitions;
  delete[] quantizedPoints;

  return 0;
}

  
