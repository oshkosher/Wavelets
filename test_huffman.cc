#include "huffman.h"


int main(int argc, char **argv) {
  
  Huffman huff(8);
  huff.update(0, 5);
  huff.update(1, 2);
  huff.update(2, 10);
  huff.update(3, 12);
  huff.update(4, 2);
  /*
  huff.update(5, 3);
  huff.update(6, 1);
  huff.update(7, 9);
  */

  huff.computeHuffmanCoding();

  huff.printEncoding();

  FILE *f = fopen("test_huffman.out", "wb");
  BitStreamFileSink fileSink(f);
  BitStreamWriter<BitStreamFileSink> bitWriter(&fileSink);
  int values[] = {0, 2, 3, 4, 2, 3, 1, 1, 1};
  int valueCount = sizeof values / sizeof(int);
  huff.encodeToStream(&bitWriter, values, valueCount);
  bitWriter.flush();
  fclose(f);

  f = fopen("test_huffman.out", "rb");
  BitStreamFileSource fileSource(f);
  BitStreamReader<BitStreamFileSource> bitReader(&fileSource);
  int *values2 = new int[valueCount];
  int readCount = huff.decodeFromStream(values2, valueCount, &bitReader);
  fclose(f);

  printf("Wrote %d, read %d\n", valueCount, readCount);
  for (int i=0; i < valueCount; i++) {
    printf("%d\t%d\n", values[i], values2[i]);
    if (values[i] != values2[i]) printf("ERROR\n");
  }
  delete[] values2;

  return 0;
}
