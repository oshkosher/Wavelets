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

static const float CDF97_SYNTHESIS_LOWPASS_FILTER[] = {
  .788485616405660,
  .377402855612650,
 -.040689417609558,
 -.02384946501938
};

static const float CDF97_SYNTHESIS_HIGHPASS_FILTER[] = {
 -.85269867900940,
  .418092273222210,
  .110624404418420,
 -.064538882628938,
 -.037828455506995
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


// If the length is not a multiple of 2^steps, it is not padded
bool is_padded_for_wavelet(int length, int steps) {
  int multiple = 1 << steps;
  return (length & (multiple - 1)) == 0;
}

bool is_padded_for_wavelet(scu_wavelet::int3 size, scu_wavelet::int3 steps) {
  return
    is_padded_for_wavelet(size.x, steps.x) &&
    is_padded_for_wavelet(size.y, steps.y) &&
    is_padded_for_wavelet(size.z, steps.z);
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

  void visitRow(NUM *row, int len, int y, int z) {
    haar_internal(len, row, steps, temp);
  }
};


template<class NUM>
void haar_3d_one_axis(CubeNum<NUM> *data, int stepCount) {
  if (stepCount > 0) {
    HaarRowVisitor<NUM> rowIter(stepCount, data->width());
    // data->template visitRows<HaarRowVisitor<NUM>>(rowIter);
    data->visitRows(rowIter);
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

  void visitRow(NUM *row, int len, int y, int z) {
    haar_inv(len, row, steps, temp);
  }
};


template<class NUM>
void haar_inv_3d_one_axis(CubeNum<NUM> *data, int stepCount) {
  if (stepCount > 0) {
    HaarInvRowVisitor<NUM> rowIter(stepCount, data->width());
    // data->template visitRows<HaarInvRowVisitor<NUM>>(rowIter);
    data->visitRows(rowIter);
  }
}



void haar_3d(CubeFloat *data, int3 stepCount, bool inverse,
             bool standardTranspose) {

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
  int length;  // length of the actual data
  const float *array;

public:
  MirroredArray(int length_, const float *array_)
    : length(length_), array(array_) {}

  float operator[] (int offset) const {

    // negative offset: mirror to a positive
    if (offset < 0) offset = -offset;

    // past the end: fold it back, repeat if necessary
    // try using modulo, see if it speeds this up
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
        

/**
   To compare with Matlab implementation in:
     Wavelets/Octave/Sergio_Matlab/fwt_1d.m

   v=[7 2 3 4 5 3 1 5 7 1 9 2 4 8 2 6]';
   udat = [v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v];
   [qa,qs] = qmfilter(' VS_7.9');
   w=fwt_1d(udat,qa,qs,1)';
   w(1,1:16)

   With a 9-element filter (4 before + center + 4 after)
     [0]: f0*x0 + 2*f1*x1 + 2*f2*x2 + 2*f3*x3 + 2*f4*x4
     [1]: f0*x1 + f1*(x0+x2) + f2*(x1+x3) + f3*(x2+x4) + f4*(x3+x5)
     [2]: f0*x2 + f1*(x1+x3) + f2*(x0+x4) + f3*(x1+x5) + f4*(x2+x6)
     [3]: f0*x3 + f1*(x2+x4) + f2*(x1+x5) + f3*(x0+x6) + f4*(x1+x7)
     [4]: f0*x4 + f1*(x3+x5) + f2*(x2+x6) + f3*(x1+x7) + f4*(x0+x8)

     [5..len-6]: normal

     [len-5]: similar
     [len-4]: 
     [len-3]: 
     [len-2]: 
     [len-1]: f0*x[n-1] + 2*f1*x[n-2] + 2*f2*x[n-3] + 2*f3*x[n-4]
                        + 2*f2*x[n-5]


    After one transform: 7.00224, 3.52709, 6.82547, 2.87912, 7.34827, 7.39307, 6.14146, 6.57857, 2.33178, -0.122068, -0.136087, -1.33848, 5.86312, 3.64358, -4.18374, -2.92383
*/
void cdf97(int length, float *data, int stepCount, float *tempGiven) {

  // check that the given array is sufficiently padded
  assert(length == dwt_padded_length(length, stepCount, false));

  float *temp;
  if (tempGiven) {
    temp = tempGiven;
  } else {
    temp = new float[length];
  }

  // encapsulate the array, automatically mirror indices past the ends
  MirroredArray array(length, data);

  int stepLength = length;
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
      temp[lowIdx++] = sum;

      arrayIdx++;

      sum = CDF97_ANALYSIS_HIGHPASS_FILTER[0] * array[arrayIdx];
      for (int filterIdx = 1; filterIdx <= 3; filterIdx++) {
        sum += CDF97_ANALYSIS_HIGHPASS_FILTER[filterIdx]
          * (array[arrayIdx + filterIdx] + array[arrayIdx - filterIdx]);
      }
      temp[hiIdx++] = sum;

      arrayIdx++;
    }

    // overwrite outputData with the results
    memcpy(data, temp, sizeof(float) * stepLength);

    /*
    printf("After step %d\n", (stepNo+1));
    for (int i=0; i < length; i++)
      printf("%f, ", (double) data[i]);
    putchar('\n');
    */
  }

  if (!tempGiven) delete[] temp;
}

/*
  Start with 7.00224, 3.52709, 6.82547, 2.87912, 7.34827, 7.39307, 6.14146, 6.57857,     2.33178, -0.122068, -0.136087, -1.33848, 5.86312, 3.64358, -4.18374, -2.92383

  Interleave: 7.00266457, 2.33238411, 3.52769017, -0.121900544, 6.82611752, -0.135406822, 2.87938881, -1.33843088, 7.34900379, 5.8633132, 7.3933816, 3.64345884, 6.14212084, -4.18343306, 6.57914591, -2.92348981
  e=[d[0], d[8], d[1], d[9], d[2], d[10], d[3], d[11], d[4], d[12], d[5], d[13], d[6], d[14], d[7], d[15]]

  evens = x[i]*.788 + (x[i-1] + x[i+1]) * .377 + (x[i-2] + x[i+2]) * -0.040 + (x[i-3] + x[i+3]) * -0.023

  odds = same, but -0.852, 0.418, 0.110, -0.0645, -0.0378
*/
void cdf97_inverse(int length, float *data, int stepCount, float *tempGiven) {
  assert(length == dwt_padded_length(length, stepCount, false));

  float *temp;
  if (tempGiven) {
    temp = tempGiven;
  } else {
    temp = new float[length];
  }

  // encapsulate the array, automatically mirror indices past the ends
  MirroredArray array(length, temp);
  
  int stepLength = length >> (stepCount-1);
  for (int stepNo=0; stepNo < stepCount; stepNo++, stepLength *= 2) {
    
    int half = stepLength >> 1;
    array.setLength(stepLength);

    // interleave: 01234567 -> 04152636
    for (int j=0; j < half; j++) {
      temp[j*2] = data[j];
      temp[j*2+1] = data[half+j];
    }

    int i=0;
    while (i < stepLength) {
      // evens
      data[i] = array[i] * CDF97_SYNTHESIS_LOWPASS_FILTER[0];
      for (int j=1; j <= 3; j++)
        data[i] += (array[i+j] + array[i-j]) * CDF97_SYNTHESIS_LOWPASS_FILTER[j];
      i++;

      // odds
      data[i] = array[i] * CDF97_SYNTHESIS_HIGHPASS_FILTER[0];
      for (int j=1; j <= 4; j++)
        data[i] += (array[i+j] + array[i-j]) * CDF97_SYNTHESIS_HIGHPASS_FILTER[j];
      i++;
    }
  }

  if (!tempGiven) delete[] temp;
}


// This will be called once for each row in the dataset by
// Cube<NUM>::visitRows. It will perform run cdf97() on the row.
template<class NUM>
class CDF97RowVisitor {
public:
  int steps;
  NUM *temp;

  CDF97RowVisitor(int s, int rowLength) : steps(s) {
    temp = new NUM[rowLength];
  }

  ~CDF97RowVisitor() {
    delete[] temp;
  }

  void visitRow(NUM *row, int len, int y, int z) {
    cdf97(len, row, steps, temp);
  }
};


// Use Cube<NUM>::visitRows to run cdf97() on every row
template<class NUM>
void cdf97_3d_one_axis(CubeNum<NUM> *data, int stepCount) {
  if (stepCount > 0) {
    CDF97RowVisitor<NUM> rowIter(stepCount, data->width());
    data->template visitRows<CDF97RowVisitor<NUM>>(rowIter);
  }
}


// This will be called once for each row in the dataset by
// Cube<NUM>::visitRows. It will perform run cdf97_inverse() on the row.
template<class NUM>
class CDF97InvRowVisitor {
public:
  int steps;
  NUM *temp;

  CDF97InvRowVisitor(int s, int rowLength) : steps(s) {
    temp = new NUM[rowLength];
  }

  ~CDF97InvRowVisitor() {
    delete[] temp;
  }

  void visitRow(NUM *row, int len, int y, int z) {
    cdf97_inverse(len, row, steps, temp);
  }
};


// Use Cube<NUM>::visitRows to run cdf97_inverse() on every row
template<class NUM>
void cdf97_inv_3d_one_axis(CubeNum<NUM> *data, int stepCount) {
  if (stepCount > 0) {
    CDF97InvRowVisitor<NUM> rowIter(stepCount, data->width());
    data->template visitRows<CDF97InvRowVisitor<NUM>>(rowIter);
  }
}


// 3-d CDF 9.7
void cdf97_3d(CubeFloat *data, scu_wavelet::int3 stepCount, bool inverse,
              bool standardTranspose, bool quiet) {

  double xform1=0, xform2=0;

  if (!inverse) {
    
    if (standardTranspose) {

      // forward standard: x steps, transpose, y steps, transpose, z steps

      cdf97_3d_one_axis(data, stepCount.x);

      double startXform = NixTimer::time();
      data->transpose3dFwd();  // xyz -> yzx
      xform1 = NixTimer::time() - startXform;

      cdf97_3d_one_axis(data, stepCount.y);

      startXform = NixTimer::time();
      data->transpose3dFwd();  // yzx -> zxy
      xform2 = NixTimer::time() - startXform;

      cdf97_3d_one_axis(data, stepCount.z);

    } else {

      // forward nonstandard: x step, transpose, y step, transpose, z step,
      // extract (1/2)^3 partial cube
      //   x, transpose, y, transpose, z, transpose, paste into full
      // extract (1/4)^3 partial cube
      //   ...

      fprintf(stderr, "nonstandard 3d CDF97 transform not implemented yet\n");
    
    }

  } else {

    if (standardTranspose) {

      // backwards standard: z inverse steps, reverse transpose,
      // y inverse steps, reverse transpose, x inverse steps
      cdf97_inv_3d_one_axis(data, stepCount.z);

      double startXform = NixTimer::time();
      data->transpose3dBack();  // zxy -> yzx
      xform1 = NixTimer::time() - startXform;

      cdf97_inv_3d_one_axis(data, stepCount.y);

      startXform = NixTimer::time();
      data->transpose3dBack();  // yzx -> xyz
      xform2 = NixTimer::time() - startXform;

      cdf97_inv_3d_one_axis(data, stepCount.x);
   
    } else {

      fprintf(stderr, "nonstandard 3d CDF97 transform not implemented yet\n");

    }
  }

  // add silly check to disable without the compiler complaining about
  // unused variables
  if (!quiet && xform1 < 0) {
    printf("Transpose %.3f ms + %.3f ms = %.3f ms\n", 
           xform1*1000, xform2*1000, (xform1+xform2)*1000);
  }
}
