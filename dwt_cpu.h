#ifndef __DWT_CPU_H__
#define __DWT_CPU_H__

#include "wavelet.h"

/*
  Simple implementation of a discrete wavelet transform using the CPU.
*/

unsigned countLeadingZeros(unsigned x);
unsigned ceilLog2(unsigned x);

// Returns the maximum number of steps a DWT can take for a given input length
// Is essentially ceil(log2(length))
int dwtMaximumSteps(int length);


// Transpose a square matrix.
void transpose_square(int size, float data[]);
void transpose_square(int size, double data[]);

/*
  Transpose an upper-left square of a square matrix.

  transpose_square_submatrix(4, 2, data):  
    a b . .       a c . .
    c d . .   ->  b d . .
    . . . .       . . . .
    . . . .       . . . .

    All "." cells would be unchanged.
*/
void transpose_square_submatrix(int total_size, int submatrix_size,
                                float data[]);

// print a matrix (for debugging purposes)
void print_matrix(int width, int height, float *data);

// Haar wavelet filter on one row of data, and the inverse.
// stepCount is the number of passes over the data. Values <= 0 will
// result in ceil(log2(data)) passes.
void haar(int length, float data[], bool inverse = false,
          int stepCount = -1);
void haar(int length, double data[], bool inverse = false,
          int stepCount = -1);

// Haar wavelet transform. on a 2-d square of data
// Returns the time the operation took in milliseconds.
float haar_2d(int size, float *data,
              bool inverse = false, int stepCount = -1,
              bool standardTranspose = false);

float haar_2d(int size, double *data,
              bool inverse = false, int stepCount = -1,
              bool standardTranspose = false);

// 3-d Haar
float haar_3d(CubeFloat *data, scu_wavelet::int3 stepCount,
              bool inverse = false, bool standardTranspose = false);

typedef enum {
  ZERO_FILL,  // fill pad elements with zero
  REFLECT,    // fill with reflection: abcde -> abcdedcb
  REPEAT      // fill with copies of last value: abcde->abcdeeee
} DWTPadding;


/*
  Given an input length, return that length rounded up to a length
  compatible with 'stepCount' steps of discrete wavelet transforms.
  If powerOfTwo is true, round up to a power of two. Otherwise,
  round up to a multiple of 2^stepCount. Return the rounded up length.
*/
int dwt_padded_length(int length, int stepCount, bool powerOfTwo);


/*
  Pad an array to the given length with the given padding method.

  The output array is returned. If output is NULL, a new array will be
  allocated. If inputLen==0, then the output array will be zero-filled.
*/
float *dwt_pad(int inputLen, float input[], 
	       int outputLen, float *output,
	       DWTPadding pad);
double *dwt_pad(int inputLen, double input[], 
		int outputLen, double *output,
		DWTPadding pad);

float *dwt_pad_2d(int rows, int cols, int rowPitch, float *input,
		  int outputRows, int outputCols, int outputPitch,
		  float *output, DWTPadding pad);
double *dwt_pad_2d(int rows, int cols, int rowPitch, double *input,
		   int outputRows, int outputCols, int outputPitch,
		   double *output, DWTPadding pad);

// do a CDF 9.7 wavelet transform
void cdf97(int length, float *data, int stepCount, float *tempGiven = NULL);
void cdf97_inverse(int length, float *data, int stepCount, float *tempGiven = NULL);


#endif // __DWT_CPU_H__
