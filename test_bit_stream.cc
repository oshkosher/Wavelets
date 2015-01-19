#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cassert>
#include <cmath>
#include <vector>
#include "nixtimer.h"
#include "bit_stream.h"

using namespace std;

#define FILENAME "test_bit_stream.out"
// #define BitStreamWriter BitStreamWriterSimple
// #define BitStreamReader BitStreamReaderSimple

bool readFile(vector<unsigned char> &data, const char *filename) {
  FILE *inf = fopen(filename, "rb");
  if (!inf) return false;
  
  data.resize(0);
  int c;
  while ((c = fgetc(inf)) != EOF) {
    data.push_back((char)c);
  }
  fclose(inf);
  return true;
}


int decodeHexDigit(char c) {
  if (c <= '9') {
    return c - '0';
  } else {
    c = toupper(c) - 'A';
    return c + 10;
  }
}


// Decode the given 2-character hexadecimal string. If it is not a valid
// hex string, return -1.
int decodeHexByte(const char *s) {
  if (!isxdigit(s[0]) || !isxdigit(s[1])) return -1;
  return decodeHexDigit(s[0]) * 16 + decodeHexDigit(s[1]);
}


bool matchesFile(const char *hexString, const char *filename) {
  vector<unsigned char> data;

  if (!readFile(data, filename)) {
    fprintf(stderr, "Failed to read \"%s\".\n", filename);
    return false;
  }

  for (size_t i=0; i < data.size(); i++) {
    int hexChar = decodeHexByte(hexString + 2*i);
    if (hexChar != data[i]) {
      fprintf(stderr, "Mismatch at byte offset %llu %02x vs %02x\n",
	      (long long unsigned)i, hexChar, data[i]);
      return false;
    }
  }
  if (hexString[data.size()*2]) {
    fprintf(stderr, "Data length mismatch\n");
    return false;
  }

  return true;
}


BitStreamWriter<BitStreamFileSink> *createStream(FILE *outf) {
  BitStreamFileSink *sink = new BitStreamFileSink(outf);
  return new BitStreamWriter<BitStreamFileSink>(sink);
}

BitStreamReader<BitStreamFileSource> *createInStream(FILE *inf) {
  BitStreamFileSource *source = new BitStreamFileSource(inf);
  return new BitStreamReader<BitStreamFileSource>(source);
}
  

void testBitStream() {
  FILE *outf, *inf;
  BitStreamWriter<BitStreamFileSink> *stream;
  BitStreamReader<BitStreamFileSource> *in;

  outf = fopen(FILENAME, "wb");
  stream = createStream(outf);

  stream->write(3, 8);
  stream->write(5, 8);
  stream->write(7, 8);
  stream->write(15, 8);
  stream->flush();
  delete stream;
  fclose(outf);
  assert(matchesFile("0305070f", FILENAME));

  inf = fopen(FILENAME, "rb");
  in = createInStream(inf);
  assert(in->read(8) == 3);
  assert(in->read(8) == 5);
  assert(in->read(8) == 7);
  assert(in->read(8) == 15);
  delete in;
  fclose(inf);

  outf = fopen(FILENAME, "wb");
  stream = createStream(outf);
  stream->write(3, 8);
  stream->write(5, 8);
  stream->write(7, 8);
  stream->write(0xffffffff, 32);
  stream->flush();
  fclose(outf);
  delete stream;
  assert(matchesFile("030507ffffffff00", FILENAME));

  inf = fopen(FILENAME, "rb");
  in = createInStream(inf);
  assert(in->read(8) == 3);
  assert(in->read(8) == 5);
  assert(in->read(8) == 7);
  assert(in->read(32) == 0xffffffff);
  delete in;
  fclose(inf);

  printf("testBitStream OK\n");
}


void testWriteArray() {
  FILE *outf, *inf;
  BitStreamWriter<BitStreamFileSink> *stream;
  BitStreamReader<BitStreamFileSource> *in;
  unsigned tmp, array[2048];
  srand(42);
  for (int i=0; i < 2048; i++)
    array[i] = rand();
  
  outf = fopen(FILENAME, "wb");
  stream = createStream(outf);

  // add 5 bits, un-aligning the rest of the data
  tmp = 21;
  stream->write(&tmp, 5);

  // add 10 words and 20 bits
  stream->write(array, 32*10 + 20);

  // add 7 bits, aligning the data
  tmp = 0x7f;
  stream->write(&tmp, 7);

  // add 10 more words
  stream->write(array+11, 10*32);
  delete stream;
  fclose(outf);

  inf = fopen(FILENAME, "rb");
  in = createInStream(inf);
  assert(in->read(5) == 21);
  for (int i=0; i < 10; i++) {
    assert(in->read(32) == array[i]);
  }
  unsigned leftover = in->read(20);
  assert(leftover == (array[10] & 0xfffff));
  assert(in->read(7) == 0x7f);

  for (int i=11; i < 21; i++) {
    assert(in->read(32) == array[i]);
  }
  delete in;
  fclose(inf);

  printf("testWriteArray OK\n");
}


void testBitStreamRandom(bool usePointer = false) {
  FILE *outf, *inf;
  BitStreamWriter<BitStreamFileSink> *stream;
  BitStreamReader<BitStreamFileSource> *in;
  int count = 100000;

  srand(42);
  outf = fopen(FILENAME, "wb");
  stream = createStream(outf);
  for (int i=0; i < count; i++) {
    int len = rand() % 31 + 1;
    unsigned value = rand() & ((1 << len) - 1);
    // stream->write(value, len);
    if (usePointer)
      stream->write(&value, len);
    else
      stream->write(value, len);
  }
  stream->flush();
  fclose(outf);
  delete stream;

  srand(42);
  inf = fopen(FILENAME, "rb");
  in = createInStream(inf);
  for (int i=0; i < count; i++) {
    int len = rand() % 31 + 1;
    int value = rand() & ((1 << len) - 1);
    int returnedValue = in->read(len);
    if (value != returnedValue) {
      fprintf(stderr, "Random test failed on i=%d\n", i);
      assert(false);
    }
  }
  delete in;
  fclose(inf);
  printf("%d random tests OK\n", count);
}

/**
   return the length of the randomized bit string
   choose the lengths with a logarithmic distribution - mostly short ones,
   but some long. Make them all multiples of 4 so the hex strings are
   easy to check visually.
*/
int randomizeBits(unsigned array[100]) {
  // [1..1000000]
  int r1 = (rand() % 1000000) + 1;

  // (e^0 .. e^6.5] == (1..665] *4 = (4..2660]
  int len = (int) exp(6.5 * r1 / 1000000.0) * 4;
  int wordCount = (len + 31) / 32;
  for (int i=0; i < wordCount; i++)
    array[i] = rand();

  return len;
}


bool bitArraysMatch(const unsigned *array, const unsigned *array2, int len) {
  while (len >= 32) {
    if (*array++ != *array2++) return false;
    len -= 32;
  }

  if (len == 0) return true;

  unsigned mask = (1 << len) - 1;
  return (*array & mask) == (*array2 & mask);
}

void printBitArray(const unsigned *array_u, int len) {
  const unsigned char *array = (const unsigned char *) array_u;
  while (len >= 8) {
    printf("%02x", *array++);
    len -= 8;
  }
  if (len > 0) {
    unsigned mask = (1 << len) - 1;
    printf("%02x", *array & mask);
  }
  putchar('\n');
}
  

void testRandomArrays() {
  FILE *outf, *inf;
  BitStreamWriter<BitStreamFileSink> *stream;
  BitStreamReader<BitStreamFileSource> *in;
  int count = 100000;
  unsigned array[100], array2[100];

  srand(42);
  outf = fopen(FILENAME, "wb");
  stream = createStream(outf);
  for (int i=0; i < count; i++) {
    int len = randomizeBits(array);
    stream->write(array, len);
    // printf("[%d] len=%d\n", i, len);
    // printBitArray(array, len);
  }
  delete stream;
  fclose(outf);

  srand(42);
  inf = fopen(FILENAME, "rb");
  in = createInStream(inf);
  for (int i=0; i < count; i++) {
    int len = randomizeBits(array);
    in->read(array2, len);
    if (!bitArraysMatch(array, array2, len)) {
      fprintf(stderr, "Random test failed on i=%d\n", i);
      printBitArray(array, len);
      printBitArray(array2, len);
      assert(false);
    }
  }
  delete in;
  fclose(inf);
  printf("%d random array tests OK\n", count);
}


double timeWriter(size_t targetSize) {
  FILE *inf = fopen(FILENAME, "wb");
  BitStreamFileSink sink(inf);
  BitStreamWriter<BitStreamFileSink> stream(&sink);
  size_t bitsWritten = 0;
  
  double startSec = NixTimer::time();
  while (bitsWritten < targetSize) {
    stream.write(0x12345678, 2);
    stream.write(0x12345678, 3);
    stream.write(0x12345678, 5);
    stream.write(0x12345678, 7);
    stream.write(0x12345678, 10);
    stream.write(0x12345678, 13);
    stream.write(0x12345678, 20);
    stream.write(0x12345678, 27);
    bitsWritten += (2+3+5+7+10+13+20+27);
  }
  stream.flush();
  fclose(inf);
  return NixTimer::time() - startSec;
}


double timeReader(size_t targetSize) {
  FILE *inf = fopen(FILENAME, "rb");
  BitStreamFileSource source(inf);
  BitStreamReader<BitStreamFileSource> stream(&source);

  double startSec = NixTimer::time();
  unsigned bitSum = 0;
  for (size_t i = 0; i < targetSize; i++) {
    bitSum += stream.readBit();
  }

  fclose(inf);
  if (bitSum == 0) printf("Odd, bit sum = 0\n");
  return NixTimer::time() - startSec;
}

void testSpeed() {
  size_t size = 100000;
  double elapsed = 0;

  while (size < 2000000000) {
    elapsed = timeWriter(size);
    // printf("%llu %.3f\n", (long long unsigned)size, elapsed);
    if (elapsed > 5) break;
    size *= 2;
  }
  double mb = size / 8.0 / 1024 / 1024;
  printf("Write %.2f MB in %.3fs: %.2f MB/sec\n", mb, elapsed, mb / elapsed);

  size_t maxSize = size;
  size = 100000;
  while (size < maxSize) {
    elapsed = timeReader(size);
    if (elapsed > 5) break;
    size *= 2;
  }
  mb = size / 8.0 / 1024 / 1024;
  printf("Read %.2f MB in %.3fs: %.2f MB/sec\n", mb, elapsed, mb / elapsed);
}


int main() {
  testBitStream();
  testBitStreamRandom(true);
  testBitStreamRandom(false);
  testWriteArray();
  testRandomArrays();
  testSpeed();
  printf("OK\n");
  return 0;
}
  
