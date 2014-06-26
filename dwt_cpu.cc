#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "nixtimer.h"

#define SQRT2     1.4142135623730950488f
#define INV_SQRT2 0.70710678118654752440f
#define DO_NORMALIZE 0

static float haar_inv_not_lifting_2d(int width, int height, float *data);


void print_matrix(int width, int height, float *data) {
  for (int y=0; y < height; y++) {
    printf("%8.5f", *data++);
    for (int x=1; x < width; x++) {
      printf(" %8.5f", *data++);
    }
    printf("\n");
  }
}


void haar_not_lifting(int length, float data[]) {
  float *temp = new float[length];
  float *s, *d;
  unsigned sampleCount = (unsigned) length;

  s = temp;

  while (sampleCount > 1) {
    int half = sampleCount >> 1;

    d = s + half;

    for (int i=0; i < half; i++) {
      d[i] = data[2*i + 1] - data[2*i];
      s[i] = data[2*i] + .5f * d[i];

#if DO_NORMALIZE
      s[i] *= SQRT2;
      d[i] *= INV_SQRT2;
#endif

    }
    memcpy(data, temp, sizeof(float) * sampleCount);
    // print_matrix(length, 1, data);

    sampleCount = half;
  }
  delete temp;
}


void haar_inv_not_lifting(int length, float data[]) {
  float *temp = new float[length];
  float *s, *d;
  int sampleCount = 2;

  s = data;

  while (sampleCount <= length) {
    int half = sampleCount >> 1;

    d = s + half;

    for (int i=0; i < half; i++) {
#if DO_NORMALIZE
      s[i] *= INV_SQRT2;
      d[i] *= SQRT2;
#endif

      temp[2*i]     = s[i] - .5f * d[i];
      temp[2*i + 1] = temp[2*i] + d[i];
    }
    memcpy(data, temp, sizeof(float) * sampleCount);
    // print_matrix(length, 1, data);

    sampleCount <<= 1;
  }
  delete temp;
}


void transpose_inplace(int width, int height, float data[]) {
  if (height != width) {
    fprintf(stderr, "transpose_inplace() only works on square data\n");
    exit(1);
  }

  for (int y=1; y < height; y++) {
    for (int x=0; x < y; x++) {
      float *p1 = data + y*width + x, *p2 = data + x*width + y;
      float tmp = *p1;
      *p1 = *p2;
      *p2 = tmp;
    }
  }
}


void haar_lifting(int length, float data[]) {
  for (int skip = 1; skip < length; skip *= 2) {
    for (int i = 0; i < length; i += skip*2) {
      data[i+skip] -= data[i]; // diff
      data[i] += data[i+skip] / 2;  // signal
    }
    // printArray(data);
  }
}    


/*
  Run Haar discrete wavelet transform on the given matrix.
  Do inverse transform iff 'inverse' is true.
  Returns the number of milliseconds elapsed.
*/
float haar_not_lifting_2d(int width, int height, float *data, bool inverse) {

  if (inverse)
    return haar_inv_not_lifting_2d(width, height, data);

  double startSec = NixTimer::time();

  // transform rows
  float *p = data;
  for (int row=0; row < height; row++) {
    haar_not_lifting(width, p);
    p += width;
  }

  transpose_inplace(width, height, data);

  // transform columns
  p = data;
  for (int col=0; col < width; col++) {
    haar_not_lifting(height, p);
    p += height;
  }

  transpose_inplace(height, width, data);

  return 1000 * (NixTimer::time() - startSec);
}


static float haar_inv_not_lifting_2d(int width, int height, float *data) {

  double startSec = NixTimer::time();

  transpose_inplace(height, width, data);

  // transform columns
  float *p = data;
  for (int row=0; row < height; row++) {
    haar_inv_not_lifting(width, p);
    p += width;
  }

  transpose_inplace(width, height, data);

  // transform rows
  p = data;
  for (int col=0; col < width; col++) {
    haar_inv_not_lifting(height, p);
    p += height;
  }

  return 1000 * (NixTimer::time() - startSec);
}
