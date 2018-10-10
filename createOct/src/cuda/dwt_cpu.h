#ifndef __DWT_CPU_H__
#define __DWT_CPU_H__

/*
  Simple implementation of a discrete wavelet transform using the CPU.
*/

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



#endif // __DWT_CPU_H__
