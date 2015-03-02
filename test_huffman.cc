#include "huffman.h"


// read integers from stdin, encode and decode them
int main(int argc, char **argv) {
  FILE *inf = stdin;
  // FILE *inf = fopen("test_huffman.in", "rt");

  if (argc < 2 || argc > 3) {
    printf("\n  test_huffman <max> [<dupValue>]\n\n");
    return 1;
  }
  int maxValue, dupValue = -1;
  if (1 != sscanf(argv[1], "%d", &maxValue)) {
    printf("Invalid max value: %s\n", argv[1]);
    return 1;
  }
  if (argc > 2) {
    if (1 != sscanf(argv[2], "%d", &dupValue)) {
      printf("Invalid dup value: %s\n", argv[2]);
      return 1;
    }
  }

  // read all the integers in inf
  std::vector<int> inputData;
  int value;
  while (fscanf(inf, "%d", &value) == 1) {
    if (value < 0 || value > maxValue) {
      printf("Value out of range: %d\n", value);
      return 1;
    }
    inputData.push_back(value);
  }
  
  printf("%d values read\n", (int)inputData.size());

  Huffman huffman(maxValue+1);
  if (dupValue >= 0) huffman.setDuplicateValue(dupValue);

  for (size_t i=0; i < inputData.size(); i++) {

    if (inputData[i] == dupValue) {
      int startDup = i;
      while (i+1 < inputData.size() && inputData[i+1] == dupValue)
	i++;
      int dupLength = i - startDup+1;
      // printf("Dup string length %d from %d to %d\n", dupLength, startDup, (int)i);
      
      huffman.addDupString(dupLength);
    } else {
      huffman.update(inputData[i], 1);
    }
  }

  huffman.computeHuffmanCoding();

  int dupKey = huffman.getDuplicateKey();

  printf("Encoded length = %d\n", huffman.totalEncodedLengthBytes());
  // huffman.printEncoding();

  std::vector<unsigned> membuf;
  BitStreamMemorySink memsink(&membuf);
  BitStreamWriter<BitStreamMemorySink> bitWriter(&memsink);

  huffman.encodeToStream(&bitWriter, inputData.data(), inputData.size());
  printf("%d words, %d bytes used\n", (int)membuf.size(), (int)membuf.size()*4);


  std::vector<int> decodeTable;
  huffman.getDecoderTable(decodeTable);

  BitStreamMemorySource memsource(&membuf);
  BitStreamReader<BitStreamMemorySource> bitReader(&memsource);
  std::vector<int> outputData;
  outputData.resize(inputData.size());

  HuffmanDecoder decoder;
  decoder.init(decodeTable);
  if (dupValue >= 0)
    decoder.setDuplicate(dupKey, dupValue);
  decoder.decodeFromStream(outputData.data(), inputData.size(), &bitReader);

  int errs = 0;
  for (size_t i=0; i < inputData.size(); i++) {
    if (inputData[i] != outputData[i]) {
      printf("Mismatch at %d: input=%d, output=%d\n", (int)i,
	     inputData[i], outputData[i]);
      errs++;
    }
  }
  printf("%d values tested, %d errors\n", (int)inputData.size(), errs);

}


// run a built-in test
int main_2(int argc, char **argv) {
  
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
