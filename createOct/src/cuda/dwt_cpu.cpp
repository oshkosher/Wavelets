#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include "nixtimer.h"
#include "dwt_cpu.h"

#define SQRT2     1.4142135623730950488f
#define INV_SQRT2 0.70710678118654752440f


// This is different from the externally visible haar_not_lifting() function
// in that no sanity check is done on the stepCount argument.
template<typename T>
static void haar_not_lifting_internal
  (int length, T data[], int stepCount);

template<typename T>
static void haar_inv_not_lifting
  (int length, T data[], int stepCount);


template<typename T>
static float haar_inv_not_lifting_2d(int size, T *data, int stepCount);


int dwtMaximumSteps(int length) {
  int steps = 0;
  while (length >= 2) {
    steps++;
    length >>= 1;
  }
  return steps;
}


void print_matrix(int width, int height, float *data) {
  for (int y=0; y < height; y++) {
    printf("%8.4f", *data++);
    for (int x=1; x < width; x++) {
      printf(" %8.4f", *data++);
    }
    printf("\n");
  }
}


template<typename T>
void haar_not_lifting_t(int length, T data[],
                        bool inverse, int stepCount) {

  // check that stepCount is valid
  int maxSteps = dwtMaximumSteps(length);
  if (stepCount < 1 || stepCount > maxSteps)
    stepCount = maxSteps;

  if (inverse)
    haar_inv_not_lifting(length, data, stepCount);
  else
    haar_not_lifting_internal(length, data, stepCount);
}

void haar_not_lifting(int length, float data[],
                      bool inverse, int stepCount) {
  haar_not_lifting_t(length, data, inverse, stepCount);
}

void haar_not_lifting(int length, double data[],
                      bool inverse, int stepCount) {
  haar_not_lifting_t(length, data, inverse, stepCount);
}


template<typename T>
static void haar_not_lifting_internal(int length, T data[],
                                      int stepCount) {
                   
  T *temp = new T[length];
  T *s, *d;
  unsigned sampleCount = (unsigned) length;

  s = temp;

  while (sampleCount > 1 && stepCount > 0) {
    int half = sampleCount >> 1;
    stepCount--;

    d = s + half;

    for (int i=0; i < half; i++) {
      d[i] = (data[2*i] - data[2*i + 1]) * INV_SQRT2;
      s[i] = (data[2*i] + data[2*i + 1]) * INV_SQRT2;
    }
    memcpy(data, temp, sizeof(T) * sampleCount);
    // print_matrix(length, 1, data);

    sampleCount = half;
  }
  delete temp;
}


template<typename T>
static void haar_inv_not_lifting(int length, T data[], int stepCount) {
  T *temp = new T[length];
  T *s, *d;

  int sampleCount = length << (stepCount - 1);

  s = data;

  while (sampleCount <= length) {
    int half = sampleCount >> 1;

    d = s + half;

    for (int i=0; i < half; i++) {
      temp[2*i]   = INV_SQRT2 * (s[i] + d[i]);
      temp[2*i+1] = INV_SQRT2 * (s[i] - d[i]);
    }
    memcpy(data, temp, sizeof(T) * sampleCount);
    // print_matrix(length, 1, data);

    sampleCount <<= 1;
  }
  delete temp;
}


// Transpose a square matrix.
template<typename T>
static void transpose_square_t(int size, T data[]) {
  for (int y=1; y < size; y++) {
    for (int x=0; x < y; x++) {
      T *p1 = data + y*size + x, *p2 = data + x*size + y;
      T tmp = *p1;
      *p1 = *p2;
      *p2 = tmp;
    }
  }
}

void transpose_square(int size, float data[]) {
  transpose_square_t(size, data);
}

void transpose_square(int size, double data[]) {
  transpose_square_t(size, data);
}


template<typename T>
static void transpose_square_submatrix_t(int total_size, int submatrix_size,
                                         T data[]) {
  for (int y=1; y < submatrix_size; y++) {
    for (int x=0; x < y; x++) {
      T *p1 = data + y*total_size + x, *p2 = data + x*total_size + y;
      T tmp = *p1;
      *p1 = *p2;
      *p2 = tmp;
    }
  }

}

void transpose_square_submatrix(int total_size, int submatrix_size,
                                float data[]) {
  transpose_square_submatrix_t(total_size, submatrix_size, data);
}

void transpose_square_submatrix(int total_size, int submatrix_size,
                                double data[]) {
  transpose_square_submatrix_t(total_size, submatrix_size, data);
}


void haar_lifting(int length, float data[], int stepCount) {

  // check that stepCount is valid
  int maxSteps = dwtMaximumSteps(length);
  if (stepCount < 1 || stepCount > maxSteps)
    stepCount = maxSteps;

  int maxSkip = 1 << (stepCount-1);

  for (int skip = 1; skip <= maxSkip; skip *= 2) {
    for (int i = 0; i < length; i += skip*2) {
      data[i+skip] -= data[i];
      data[i] += data[i+skip] * .5f;
      data[i] *= SQRT2;
      data[i+skip] *= INV_SQRT2;
    }
    // printArray(data);
  }
}    

void haar_inv_lifting(int length, float data[], int stepCount) {

  // check that stepCount is valid
  int maxSteps = dwtMaximumSteps(length);
  if (stepCount < 1 || stepCount > maxSteps)
    stepCount = maxSteps;

  for (int skip = 1 << (stepCount-1); skip >= 1; skip >>= 1) {
    for (int i = 0; i < length; i += skip*2) {
      data[i+skip] *= SQRT2;
      data[i] *= INV_SQRT2;
      data[i] -= data[i+skip] / 2;  // signal
      data[i+skip] += data[i]; // diff
    }
    // printArray(data);
  }
}    


/*
  Run Haar discrete wavelet transform on the given matrix.
  Do inverse transform iff 'inverse' is true.
  Returns the number of milliseconds elapsed.
*/
template<typename T>
static float haar_not_lifting_2d_t(int size, T *data, bool inverse,
                                   int stepCount) {

  // check that stepCount is valid
  int maxSteps = dwtMaximumSteps(size);
  if (stepCount < 1 || stepCount > maxSteps)
    stepCount = maxSteps;

  if (inverse)
    return haar_inv_not_lifting_2d(size, data, stepCount);

  double startSec = NixTimer::time();

  int stepLength = size;
  for (int step = 0; step < stepCount; step++) {

    // transform rows
    T *p = data;
    for (int row=0; row < stepLength; row++) {
      haar_not_lifting_internal(stepLength, p, 1);
      p += size;
    }

    transpose_square_submatrix(size, stepLength, data);

    // transform columns
    p = data;
    for (int col=0; col < stepLength; col++) {
      haar_not_lifting_internal(stepLength, p, 1);
      p += size;
    }

    transpose_square_submatrix(size, stepLength, data);

    stepLength >>= 1;
  }

  return (float) (1000 * (NixTimer::time() - startSec));
}

float haar_not_lifting_2d(int size, float *data, bool inverse,
                          int stepCount) {
  return haar_not_lifting_2d_t(size, data, inverse, stepCount);
}

float haar_not_lifting_2d(int size, double *data, bool inverse,
                          int stepCount) {
  return haar_not_lifting_2d_t(size, data, inverse, stepCount);
}


template<typename T>
static float haar_inv_not_lifting_2d(int size, T *data,
                                     int stepCount) {

  double startSec = NixTimer::time();

  for (int step = stepCount; step >= 1; step--) {
    int stepLength = size >> (step-1);
    
    transpose_square_submatrix(size, stepLength, data);

    // transform columns
    T *p = data;
    for (int row=0; row < stepLength; row++) {
      haar_inv_not_lifting(stepLength, p, 1);
      p += size;
    }

    transpose_square_submatrix(size, stepLength, data);

    // transform rows
    p = data;
    for (int col=0; col < stepLength; col++) {
      haar_inv_not_lifting(stepLength, p, 1);
      p += size;
    }
  }

  return (float) (1000 * (NixTimer::time() - startSec));
}
