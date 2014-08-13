#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include "nixtimer.h"
#include "dwt_cpu.h"

#define SQRT2     1.4142135623730950488f
#define INV_SQRT2 0.70710678118654752440f

#define MIN(x,y) (((x) < (y)) ? (x) : (y))

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


// return the number of high-order zero bits in a 32 bit integer
unsigned countLeadingZeros(unsigned x) {
  unsigned y, n;
  n = 32;
  y = x >> 16; if (y != 0) {n-=16; x=y;}
  y = x >>  8; if (y != 0) {n-= 8; x=y;}
  y = x >>  4; if (y != 0) {n-= 4; x=y;}
  y = x >>  2; if (y != 0) {n-= 2; x=y;}
  y = x >>  1; if (y != 0) return n - 2;
  return n - x;
}


// Return ceil(log2(x)), or the log2 of the smallest power of two
// that is greater than or equal to x.
unsigned ceilLog2(unsigned x) {
  if (x == 0) return 0;
  return 32 - countLeadingZeros(x-1);
}


// Returns the maximum number of steps a DWT can take for a given input length
int dwtMaximumSteps(int length) {
  return ceilLog2(length) - 1;
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
  delete[] temp;
}


/*
  Given an input length, return that length rounded up to a length
  compatible with 'stepCount' steps of discrete wavelet transforms.
  If powerOfTwo is true, round up to a power of two. Otherwise,
  round up to a multiple of 2^stepCount. Return the rounded up length.
*/
int dwt_padded_length(int length, int stepCount, bool powerOfTwo) {
  if (powerOfTwo) {
    return 1 << ceilLog2(length);
  } else {
    int mod = 1 << stepCount;
    return mod * ((length-1) / mod + 1);
  }
}


/*
  Pad an array to the given length with the given padding method.

  The output array is returned. If output is NULL, a new array will be
  allocated. If inputLen==0, then the output array will be zero-filled.
*/
template<typename T>
static T *dwt_pad_t(int inputLen, T input[], 
		    int outputLen, T *output,
		    DWTPadding pad) {

  // allocate output array, if necessary
  if (output == NULL) {
    output = new T[outputLen];
  }

  // if the input array is not the same as the output array, copy the first
  // 'inputLen' elements to the output array.
  if (input != output) {
    for (int i=0; i < inputLen; i++) output[i] = input[i];
  }

  // If the length is zero, there's nothing to copy, so might as well
  // fill it with zeros. If there's only one element, reflecting is the
  // same as repeating, so just repeat.
  if (inputLen == 0) {
    pad = ZERO_FILL;
  } else if (inputLen == 1 && pad == REFLECT) {
    pad = REPEAT;
  }

  switch (pad) {
  case ZERO_FILL:
    for (int i=inputLen; i < outputLen; i++) {
      output[i] = 0;
    }
    break;

  case REPEAT:
    T repeatValue;
    repeatValue = output[inputLen-1];
    for (int i=inputLen; i < outputLen; i++) {
      output[i] = repeatValue;
    }
    break;

  case REFLECT:
    // if the input array is less than half the length of the output array,
    // multiple passes will be needed. For example:
    // input: abc, output length 8
    // pass 1: abcba
    // pass 2: abcbabcb
    int writePos = inputLen;
    do {
      int iters = MIN(inputLen-1, outputLen - inputLen);
      int readPos = inputLen-2;
      for (int i=0; i < iters; i++) {
	output[writePos++] = output[readPos--];
      }
      inputLen = writePos;
    } while (writePos < outputLen);
    break;
  }
    

  return output;
}


float *dwt_pad(int inputLen, float input[], 
	       int outputLen, float *output,
	       DWTPadding pad) {
  return dwt_pad_t(inputLen, input, outputLen, output, pad);
}

double *dwt_pad(int inputLen, double input[], 
		int outputLen, double *output,
		DWTPadding pad) {
  return dwt_pad_t(inputLen, input, outputLen, output, pad);
}


template<typename T>
static void copyRow(T *dest, int destRow, int destPitch,
		    T *src,  int srcRow,  int srcPitch,
		    int cols) {
  memcpy(dest + destRow * destPitch,
	 src  + srcRow  * srcPitch,
	 sizeof(T) * cols);
}

template<typename T>
static T *dwt_pad_2d_t(int inputRows,  int inputCols, 
		     int inputPitch, T *input,
		     int outputRows, int outputCols,
		     int outputPitch, T *output, 
		     DWTPadding pad) {

  // allocate output array, if necessary
  if (output == NULL) {
    output = new T[outputRows * outputPitch];
  }

  // if the input array is not the same as the output array, copy the first
  // 'inputLen' elements to the output array.
  if (input != output) {
    for (int i=0; i < inputRows; i++) {
      copyRow(output, i, outputPitch, input, i, inputPitch, inputCols);
    }
  }
  
  // pad each row
  for (int i=0; i < inputRows; i++) {
    dwt_pad(inputCols,  input  + i*inputPitch,
	    outputCols, output + i*outputPitch, pad);
  }
  
  // pad columns by copying rows

  // If the length is zero, there's nothing to copy, so might as well
  // fill it with zeros. If there's only one element, reflecting is the
  // same as repeating, so just repeat.
  if (inputRows == 0) {
    pad = ZERO_FILL;
  } else if (inputRows == 1 && pad == REFLECT) {
    pad = REPEAT;
  }

  switch (pad) {
  case ZERO_FILL:
    for (int i=inputRows; i < outputRows; i++) {
      // fill the row with zeros
      memset(output + i*outputPitch, 0, sizeof(T) * outputCols);
    }
    break;

  case REPEAT:
    for (int i=inputRows; i < outputRows; i++) {
      copyRow(output, i, outputPitch,
	      output, inputRows-1, outputPitch, outputCols);
    }
    break;

  case REFLECT:
    // if the input array is less than half the length of the output array,
    // multiple passes will be needed. For example:
    // input: abc, output length 8
    // pass 1: abcba
    // pass 2: abcbabcb
    int inputLen = inputRows;
    int writePos = inputLen;
    do {
      int iters = MIN(inputLen-1, outputRows - inputLen);
      int readPos = inputLen-2;
      for (int i=0; i < iters; i++) {
	copyRow(output, writePos++, outputPitch,
		output, readPos--,  outputPitch, outputCols);
      }
      inputLen = writePos;
    } while (writePos < outputRows);
    break;
  }

  return output;
}

float *dwt_pad_2d(int inputRows,  int inputCols, 
		     int inputPitch, float *input,
		     int outputRows, int outputCols,
		     int outputPitch, float *output, 
		     DWTPadding pad) {
  return dwt_pad_2d_t(inputRows, inputCols, inputPitch, input,
		      outputRows, outputCols, outputPitch, output, pad);
}

double *dwt_pad_2d(int inputRows,  int inputCols, 
		     int inputPitch, double *input,
		     int outputRows, int outputCols,
		     int outputPitch, double *output, 
		     DWTPadding pad) {
  return dwt_pad_2d_t(inputRows, inputCols, inputPitch, input,
		      outputRows, outputCols, outputPitch, output, pad);
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


