#include <cstdio>
#include "lloyds.h"


int main() {
  
  float points[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  int psize = sizeof points / sizeof(float);
  float codebook[] = {-5, -1, 1, 5};
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

  
