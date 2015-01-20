#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include "nixtimer.h"
#include "dwt_cpu.h"

using namespace scu_wavelet;

#define SQRT2     1.4142135623730950488f
#define INV_SQRT2 0.70710678118654752440f

#define MIN(x,y) (((x) < (y)) ? (x) : (y))

static const float CDF97_ANALYSIS_LOWPASS_FILTER[] = {
   .85269867900940,
   .377402855612650,
  -.110624404418420,
  -.02384946501938,
   .037828455506995
};

static const float CDF97_ANALYSIS_HIGHPASS_FILTER[] = {
  -.788485616405660,
   .418092273222210,
   .040689417609558,
  -.064538882628938
};


// This is different from the externally visible haar() function
// in that no sanity check is done on the stepCount argument.
// temp is optional - if the caller wants to speed up the computation,
// it can supply an array of length 'length' so these won't need to
// allocate one.
template<typename T>
static void haar_internal
  (int length, T data[], int stepCount, T *temp = NULL);

template<typename T>
static void haar_inv
(int length, T data[], int stepCount, T *temp = NULL);


template<typename T>
static float haar_inv_2d(int size, T *data, int stepCount, bool standard);


// void haar_3d_one_axis(CubeFloat *data, int stepCount);
// void haar_inv_3d_one_axis(CubeFloat *data, int stepCount);



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
  return ceilLog2(length);
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
void haar_tmpl(int length, T data[],
            bool inverse, int stepCount) {

  // check that stepCount is valid
  int maxSteps = dwtMaximumSteps(length);
  if (stepCount < 1 || stepCount > maxSteps)
    stepCount = maxSteps;

  if (inverse)
    haar_inv(length, data, stepCount);
  else
    haar_internal(length, data, stepCount);
}

void haar(int length, float data[],
          bool inverse, int stepCount) {
  haar_tmpl(length, data, inverse, stepCount);
}

void haar(int length, double data[],
          bool inverse, int stepCount) {
  haar_tmpl(length, data, inverse, stepCount);
}


template<typename T>
static void haar_internal(int length, T data[], int stepCount, T *tempGiven) {

  T *temp;
  if (tempGiven) {
    temp = tempGiven;
  } else {
    temp = new T[length];
  }

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

  if (!tempGiven) delete[] temp;
}


template<typename T>
static void haar_inv(int length, T data[], int stepCount, T *tempGiven) {
  T *temp;
  if (tempGiven) {
    temp = tempGiven;
  } else {
    temp = new T[length];
  }
  T *s, *d;

  int sampleCount = length >> (stepCount - 1);

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
  if (!tempGiven) delete[] temp;
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


// Transpose a square matrix.
template<typename T>
static void transpose_square_tmpl(int size, T data[]) {
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
  transpose_square_tmpl(size, data);
}

void transpose_square(int size, double data[]) {
  transpose_square_tmpl(size, data);
}


template<typename T>
static void transpose_square_submatrix_tmpl(int total_size, int submatrix_size,
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
  transpose_square_submatrix_tmpl(total_size, submatrix_size, data);
}

void transpose_square_submatrix(int total_size, int submatrix_size,
                                double data[]) {
  transpose_square_submatrix_tmpl(total_size, submatrix_size, data);
}


/*
  Run Haar discrete wavelet transform on the given matrix.
  Do inverse transform iff 'inverse' is true.
  Returns the number of milliseconds elapsed.
*/
template<typename T>
static float haar_2d_tmpl(int size, T *data, bool inverse, int stepCount,
                          bool standardTranspose = false) {

  // check that stepCount is valid
  int maxSteps = dwtMaximumSteps(size);
  if (stepCount < 1 || stepCount > maxSteps)
    stepCount = maxSteps;

  if (inverse)
    return haar_inv_2d(size, data, stepCount, standardTranspose);

  double startSec = NixTimer::time();

  if (standardTranspose) {

    for (int i = 0; i < 2; i++) {

      int stepLength = size;
      for (int step = 0; step < stepCount; step++) {

        // transform rows
        T *p = data;
        for (int row=0; row < size; row++) {
          haar_internal(stepLength, p, 1);
          p += size;
        }

        stepLength >>= 1;
      }

      transpose_square(size, data);
    }

  } else {

    int stepLength = size;
    for (int step = 0; step < stepCount; step++) {

      // transform rows
      T *p = data;
      for (int row=0; row < stepLength; row++) {
        haar_internal(stepLength, p, 1);
        p += size;
      }

      transpose_square_submatrix(size, stepLength, data);

      // transform columns
      p = data;
      for (int col=0; col < stepLength; col++) {
        haar_internal(stepLength, p, 1);
        p += size;
      }

      transpose_square_submatrix(size, stepLength, data);

      stepLength >>= 1;
    }

  }

  return (float) (1000 * (NixTimer::time() - startSec));
}

float haar_2d(int size, float *data, bool inverse,
              int stepCount, bool standardTranspose) {
  return haar_2d_tmpl(size, data, inverse, stepCount, standardTranspose);
}

float haar_2d(int size, double *data, bool inverse,
              int stepCount, bool standardTranspose) {
  return haar_2d_tmpl(size, data, inverse, stepCount, standardTranspose);
}


template<typename T>
static float haar_inv_2d(int size, T *data,
                         int stepCount, bool standard) {

  double startSec = NixTimer::time();

  if (standard) {

    for (int i=0; i < 2; i++) {

      transpose_square(size, data);

      for (int step = stepCount; step >= 1; step--) {
        int stepLength = size >> (step-1);

        // transform columns
        T *p = data;
        for (int row=0; row < size; row++) {
          haar_inv(stepLength, p, 1);
          p += size;
        }

      }
    }

  } else {

    for (int step = stepCount; step >= 1; step--) {
      int stepLength = size >> (step-1);
    
      transpose_square_submatrix(size, stepLength, data);

      // transform columns
      T *p = data;
      for (int row=0; row < stepLength; row++) {
        haar_inv(stepLength, p, 1);
        p += size;
      }

      transpose_square_submatrix(size, stepLength, data);

      // transform rows
      p = data;
      for (int col=0; col < stepLength; col++) {
        haar_inv(stepLength, p, 1);
        p += size;
      }
    }
  }

  return (float) (1000 * (NixTimer::time() - startSec));
}


template<class NUM>
class HaarRowVisitor {
public:
  int steps;
  NUM *temp;

  HaarRowVisitor(int s, int rowLength) : steps(s) {
    temp = new NUM[rowLength];
  }

  ~HaarRowVisitor() {
    delete[] temp;
  }

  void visitRow(NUM *row, int len) {
    haar_internal(len, row, steps, temp);
  }
};


template<class NUM>
void haar_3d_one_axis(CubeNum<NUM> *data, int stepCount) {
  if (stepCount > 0) {
    HaarRowVisitor<NUM> rowIter(stepCount, data->width());
    data->template visitRows<HaarRowVisitor<NUM>>(rowIter);
  }
}


template<class NUM>
class HaarInvRowVisitor {
public:
  int steps;
  NUM *temp;

  HaarInvRowVisitor(int s, int rowLength) : steps(s) {
    temp = new NUM[rowLength];
  }

  ~HaarInvRowVisitor() {
    delete[] temp;
  }

  void visitRow(NUM *row, int len) {
    haar_inv(len, row, steps, temp);
  }
};


template<class NUM>
void haar_inv_3d_one_axis(CubeNum<NUM> *data, int stepCount) {
  if (stepCount > 0) {
    HaarInvRowVisitor<NUM> rowIter(stepCount, data->width());
    data->template visitRows<HaarInvRowVisitor<NUM>>(rowIter);
  }
}



float haar_3d(CubeFloat *data, int3 stepCount, bool inverse,
              bool standardTranspose) {
  double startTime = NixTimer::time();

  if (!inverse) {
    
    if (standardTranspose) {

      // forward standard: x steps, transpose, y steps, transpose, z steps

      haar_3d_one_axis(data, stepCount.x);
      data->transpose3dFwd();  // xyz -> yzx
      haar_3d_one_axis(data, stepCount.y);
      data->transpose3dFwd();
      haar_3d_one_axis(data, stepCount.z);

    } else {

      // forward nonstandard: x step, transpose, y step, transpose, z step,
      // extract (1/2)^3 partial cube
      //   x, transpose, y, transpose, z, transpose, paste into full
      // extract (1/4)^3 partial cube
      //   ...

      fprintf(stderr, "nonstandard 3d Haar transform not implemented yet\n");
    
    }

  } else {

    if (standardTranspose) {

      // backwards standard: z inverse steps, reverse transpose,
      // y inverse steps, reverse transpose, x inverse steps
      haar_inv_3d_one_axis(data, stepCount.z);
      data->transpose3dBack();
      haar_inv_3d_one_axis(data, stepCount.y);
      data->transpose3dBack();
      haar_inv_3d_one_axis(data, stepCount.x);
   
    } else {

      fprintf(stderr, "nonstandard 3d Haar transform not implemented yet\n");

    }
  }

  return (NixTimer::time() - startTime) * 1000;
}




/*
  Wrap an array such that if you reference values beyond the ends
  of the array, the results will be mirrored array values.

  For example, given an array with 7 elements: 0 1 2 3 4 5 6
  Request array[0..6] and you'll get the usual values array[0..6].
  array[-1] returns array[1], array[-2] return array[2], etc.
  array[7] returns array[5], array[8] return array[4], etc.
*/
class MirroredArray {
  float *array;
  int length;  // length of the actual data

public:
  MirroredArray(float *array_, int length_)
    : array(array_), length(length_) {}

  float operator[] (int offset) const {

    // negative offset: mirror to a positive
    if (offset < 0) offset = -offset;

    // past the end: fold it back, repeat if necessary
    while (offset >= length) {
      offset = length*2 - offset - 2;

      if (offset < 0) offset = -offset;
    }

    return array[offset];
  }

  void setLength(int len) {length = len;}
};


// Like MirroredArray, but simpler. Ask for an invalid index and you get 0.
class ZeroExtendedArray {
  float *array;
  int length;  // length of the actual data

public:
  ZeroExtendedArray(float *array_, int length_)
    : array(array_), length(length_) {}

  float operator[] (int offset) const {

    if (offset < 0 || offset >= length) return 0;

    return array[offset];
  }

  void setLength(int len) {length = len;}
};
        

void cdf97(int inputLength, const float inputData[], int stepCount,
           int *resultLength, float **resultData) {

  int len = dwt_padded_length(inputLength, stepCount, false);
  float *outputData, *tmpData = new float[len];

  if (*resultData) {
    outputData = *resultData;
  } else {
    *resultData = outputData = new float[len];
  }

  // how many cells are inserted before the original data
  int inputDataOffset = (len - inputLength) / 2;

  // copy input data into place
  memcpy(outputData+inputDataOffset, inputData, sizeof(float) * len);

  // replicate the first element
  for (int i=0; i < inputDataOffset; i++)
    outputData[i] = inputData[0];

  // replicate the last element
  for (int i=inputLength+inputDataOffset; i < len; i++)
    outputData[i] = inputData[inputLength-1];

  // encapsulate the array, automatically mirror indices past the ends
  MirroredArray array(outputData, len);
  // ZeroExtendedArray array(outputData, len);

  int stepLength = len;
  for (int stepNo=0; stepNo < stepCount; stepNo++, stepLength /= 2) {

    array.setLength(stepLength);

    int arrayIdx = 0, lowIdx = 0, hiIdx = stepLength/2;
    while (arrayIdx < stepLength) {

      // Apply the low pass filter convolution.
      // It's symmetric, so apply coefficient N to array[i+N] and array[i-N].
      float sum = CDF97_ANALYSIS_LOWPASS_FILTER[0] * array[arrayIdx];
      for (int filterIdx = 1; filterIdx <= 4; filterIdx++) {
        sum += CDF97_ANALYSIS_LOWPASS_FILTER[filterIdx]
          * (array[arrayIdx + filterIdx] + array[arrayIdx - filterIdx]);
      }
      tmpData[lowIdx++] = sum;

      arrayIdx++;

      sum = CDF97_ANALYSIS_HIGHPASS_FILTER[0] * array[arrayIdx];
      for (int filterIdx = 1; filterIdx <= 3; filterIdx++) {
        sum += CDF97_ANALYSIS_HIGHPASS_FILTER[filterIdx]
          * (array[arrayIdx + filterIdx] + array[arrayIdx - filterIdx]);
      }
      tmpData[hiIdx++] = sum;

      arrayIdx++;
    }

    // overwrite outputData with the results
    memcpy(outputData, tmpData, sizeof(float) * stepLength);

    printf("After step %d\n", stepNo);
    for (int i=0; i < len; i++)
      printf("%f, ", outputData[i]);
    putchar('\n');
  }

  delete[] tmpData;

  *resultLength = len;
  *resultData = outputData;
}

