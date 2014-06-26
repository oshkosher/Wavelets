#ifndef __DWT_CPU_H__
#define __DWT_CPU_H__

/*
  Simple implementation of a discrete wavelet transform using the CPU.
*/

// transpose a square matrix in place
void transpose_inplace(int width, int height, float data[]);

// print a matrix (for debugging purposes)
void print_matrix(int width, int height, float *data);

// Haar wavelet filter on one row of data, and the inverse.
void haar_not_lifting(int length, float data[]);
void haar_inv_not_lifting(int length, float data[]);

// Haar wavelet filter on multiple rows of data
float haar_not_lifting_2d(int width, int height, float *data,
                          bool inverse = false);

// Lifting implementation (not tested much)
void haar_lifting(int length, float data[]);




#endif // __DWT_CPU_H__
