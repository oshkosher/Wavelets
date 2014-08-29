#include <cstdio>
#include <cstdlib>
#include "data_io.h"
#include "rle.h"
#include "huffman.h"
// #include "rle_virtual.h"

/* Using an unordered_map (hash): 4.2s
   Using a map (tree): 6.2s
   Using an array: 1.4s
*/


void printHelp() {
  printf("\n  test_rle <data file>\n\n");
  exit(1);
}

// this takes about .57 seconds to run and just scans the data
// the version using virtual functions takes about 1.8 seconds.
// the version using templates uses about 1.4 seconds
int main_dummy(int argc, char **argv) {
  if (argc != 2) printHelp();

  char *inputFile = argv[1];
  float *data;
  int width, height;
  if (!readDataFile(inputFile, &data, &width, &height)) return 1;

  float sum = 0, *end = data + width*height;
  for (float *p = data; p < end; p++)
    sum += *p;

  printf("sum: %f\n", sum);

  return 0;
}

int main(int argc, char **argv) {
  if (argc != 2) printHelp();

  char *inputFile = argv[1];
  float *data;
  int width, height;
  if (!readDataFile(inputFile, &data, &width, &height)) return 1;

  //PrintValues printValues;

  /*
  FrequencyCountPair frequencyCount;
  EncodeRunLength encoder(frequencyCount);
  FloatMatrixToGrayscale matrixScan(width, height, data, encoder);
  */

  // FrequencyCountPair frequencyCount;
  // EncodeRunLength<FrequencyCountPair> encoder;
  FloatMatrixToGrayscale<EncodeRunLength<FrequencyCountPair> >
    matrixScan(width, height, data);

  matrixScan.scan();

  /*
  FrequencyCountPair::MapType::iterator iter;
  for (iter = frequencyCount.counts.begin();
       iter != frequencyCount.counts.end();
       iter++) {
    printf("%d: %d\n", iter->first, iter->second);
  }
  */

  FrequencyCountPair &frequencyCount = matrixScan.sink.sink;
  Huffman huff(256);

  for (int value=0; value < 256; value++) {
    int count = frequencyCount.counts[value];
    // printf("%d, %d\n", i, count);
    huff.update(value, count);
  }

  huff.computeHuffmanCoding();

  long normalBits = 0, encodedBits = 0;
  for (int i=0; i < 256; i++) {
    // printf("%3d: %d (%d)\n", i, huff.encodedLength(i), huff.getCount(i));
    normalBits += 8*huff.getCount(i);
    encodedBits += huff.encodedLength(i)*huff.getCount(i);
  }
  printf("flat: %ld\n", normalBits);
  printf("huff: %ld\n", encodedBits);

  delete[] data;

  return 0;
}
