#ifndef __DATA_IO__
#define __DATA_IO__

/*
  Read and write 2D arrays of floats.

  Two formats are supported: text and binary.
  When reading a file, the code will automatically detect which one
  the data is in. When writing, the user can specify which format to
  use with the "isBinary" flag.

  Text format:

    First row: <C> cols <R> rows 
    After that, <R> lines of data, where each line is one row of data
    containing <C> floating point values separated by spaces.

  Binary format:

    First 8 bytes: "binary \n"
      This is used to distinguish between the text and binary formats.
    Next 8 bytes: two 4-byte little-endian integers, <width> and <height>
    Then <width>*<height> little-endian 4-byte floats in row-major order.
*/

// Read a 2-d array of float data in either of the formats listed above.
// On error, complain to stderr and return false. *data is allocated
// with new[]. The caller is responsible for calling delete[] on it.
bool readDataFile(const char *filename, float **data, int *width, int *height);

// Same as the function above, but returns an array of doubles.
// NOTE: the data is currenly still stored as floats.
bool readDataFile(const char *filename, double **data, 
                  int *width, int *height);

// Write a 2-d array of float data. 'isBinary' specifies whether the
// to use the binary format or not.
bool writeDataFile(const char *filename, float *data, int width, int height,
                   bool isBinary = true);

// Same as above, but takes an array of doubles.
// NOTE: the data is currenly still stored as floats.
bool writeDataFile(const char *filename, double *data, int width, int height,
                   bool isBinary = true);
  


#endif // __DATA_IO__
