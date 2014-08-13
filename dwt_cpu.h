#ifndef __DWT_CPU_H__
#define __DWT_CPU_H__

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
void haar_not_lifting(int length, float data[], bool inverse = false,
                      int stepCount = -1);
void haar_not_lifting(int length, double data[], bool inverse = false,
                      int stepCount = -1);

// Haar wavelet filter on a 2-d square of data
float haar_not_lifting_2d(int size, float *data,
                          bool inverse = false, int stepCount = -1);

float haar_not_lifting_2d(int size, double *data,
                          bool inverse = false, int stepCount = -1);

// Lifting implementation (not tested much)
void haar_lifting(int length, float data[], int stepCount = -1);
void haar_inv_lifting(int length, float data[], int stepCount = -1);

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


#endif // __DWT_CPU_H__
