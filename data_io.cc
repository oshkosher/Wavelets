#include <cstdio>
#include <cstring>
#include <cassert>
#include <vector>
#include "data_io.h"

// See "data_io.h" for details on the file format.

// check if the given file is in binary format
static bool isBinary(const char *filename) {
  FILE *inf = fopen(filename, "rb");
  if (!inf) return false;

  char buf[9];
  fread(buf, 1, 8, inf);
  fclose(inf);
  buf[8] = 0;
  return !strncmp(buf, "binary", 6);
}

// read a file in binary format
static bool readBinaryDataFile(const char *filename, float **data,
                               int *width, int *height) {

  FILE *inf = fopen(filename, "rb");
  if (!inf) {
    fprintf(stderr, "Failed to open \"%s\" for reading.\n", filename);
    return false;
  }

  // skip the header
  fseek(inf, 8, SEEK_SET);
  fread(width, sizeof(int), 1, inf);
  fread(height, sizeof(int), 1, inf);

  // check for invalid image sizes
  if (*width < 1 || *height < 1) {
    fprintf(stderr, "Invalid image size: %d x %d\n", *width, *height);
    fclose(inf);
    return false;
  }

  long long unsigned cells = *height * *width;
  if (cells > 0xffffffff) {
    fprintf(stderr, "Image size too large: %d x %d\n", *width, *height);
    fclose(inf);
    return false;
  }

  *data = new float[(int)cells];
  if (!*data) {
    fprintf(stderr, "Out of memory allocating %d x %d image\n",
            *width, *height);
    fclose(inf);
    return false;
  }

  fread(*data, 4, (int)cells, inf);
  fclose(inf);

  return true;
}


static bool readTextDataFile(const char *filename, float **data,
                             int *width, int *height) {

  FILE *inf = fopen(filename, "rt");
  if (!inf) {
    fprintf(stderr, "Failed to open \"%s\" for reading.\n", filename);
    return false;
  }

  // read the first line
  if (2 != fscanf(inf, "%d cols %d rows", width, height)
      || *width < 1 || *height < 1) {
    fprintf(stderr, "Bad input file format (width=%d, height=%d)\n",
            *width, *height);
    fclose(inf);
    return false;
  }

  int entryCount = (*width) * (*height);

  *data = new float[entryCount];

  for (int i=0; i < entryCount; i++) {
    float tmp;
    if (1 != fscanf(inf, "%f", &tmp)) {
      fprintf(stderr, "Bad data: expected %d values, got %d\n", entryCount, i);
      fclose(inf);
      return false;
    }

    (*data)[i] = tmp;
  }

  fclose(inf);
  return true;
}


// Read a file of data, automatically detecting whether it's in
// text or binary format.
// On error, print it to stderr and return false.
// Caller is responsible for calling 'delete[]' on *data.
bool readDataFile(const char *filename, float **data, int *width, int *height) {
  bool success;

  if (isBinary(filename)) {
    success = readBinaryDataFile(filename, data, width, height);
  } else {
    success = readTextDataFile(filename, data, width, height);
  }

  return success;
}


// Read data as doubles.
bool readDataFile(const char *filename, double **data,
                  int *width, int *height) {
  bool success;
  float *floatData;

  if (isBinary(filename)) {
    success = readBinaryDataFile(filename, &floatData, width, height);
  } else {
    success = readTextDataFile(filename, &floatData, width, height);
  }

  if (!success) return false;

  int size = *width * *height;
  *data = new double[size];
  for (int i=0; i < size; i++)
    (*data)[i] = floatData[i];

  delete[] floatData;

  return true;
}


// Write a data file in binary format
// Note: this code assumes the current machine is little endian.
// If this is ported to a big endian machine it will need to be tweaked.
static bool writeBinaryDataFile(const char *filename, float *data, 
                                int width, int height) {

  static unsigned endianCheck = 0x11223344;
                         
  FILE *outf = fopen(filename, "wb");
  if (!outf) {
    fprintf(stderr, "Failed to open \"%s\" for writing.\n", filename);
    return false;
  }

  assert(sizeof(int) == 4);
  assert(sizeof(float) == 4);
  assert(*((char*)&endianCheck) == 0x44);

  // 8 byte header
  fwrite("binary \n", 1, 8, outf);

  // width and height
  fwrite(&width, sizeof(int), 1, outf);
  fwrite(&height, sizeof(int), 1, outf);

  size_t bytesWritten = fwrite(data, sizeof(float), width*height, outf);
  if (bytesWritten != (size_t)(width*height)) {
    fprintf(stderr, "Failed to write to \"%s\".\n", filename);
    fclose(outf);
    return false;
  }
  fclose(outf);
  return true;
}

bool writeDataFile(const char *filename, double *dataDoubles,
                   int width, int height, bool isBinary) {
                   
  float *data = new float[width*height];
  for (int i=0; i < width*height; i++)
    data[i] = (float)dataDoubles[i];
  bool success = writeDataFile(filename, data, width, height, isBinary);
  delete[] data;
  return success;
}


// Write a data file in text format
bool writeDataFile(const char *filename, float *data, int width, int height,
                   bool isBinary) {
  if (isBinary) return writeBinaryDataFile(filename, data, width, height);

  FILE *outf = fopen(filename, "wt");
  if (!outf) {
    fprintf(stderr, "Failed to open \"%s\" for writing.\n", filename);
    return false;
  }

  fprintf(outf, "%d cols %d rows\n", width, height);
  for (int r=0; r < height; r++) {
    for (int c=0; c < width; c++) {
      if (c > 0) fputc(' ', outf);
      fprintf(outf, "%7.5f", *data++);
    }
    fputc('\n', outf);
  }

  fclose(outf);
  return true;
}
