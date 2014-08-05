#include <octave/oct.h>

int haar(float *output, float *input, int width, int steps, bool inverse, int blockSize);

DEFUN_DLD (cudahaar, args, , "\n"
          "  cudahaar(input, steps = 1, inverse = false, text = false, blockSize = -1)\n"
          "  Do a Haar discrete wavelet transform.\n"
          "  Options:\n"
          "    -inverse : invert the transform\n"
          "    -text : output data in text format rather than binary\n"
          "    -blocksize : specify the thread block size\n"
          "  By default, one transformation step will be done.\n"
          "\n") {

  int nargin = args.length();

  if (nargin != 1) {
    print_usage();
  } else {
    int steps = 1;
    bool inverse = false;
    bool textOutput = false;
    int blockSize = -1;

    if (nargin >= 2) {steps  = args(1).int_value();}
    if (nargin >= 3) {inverse    = args(2).bool_value();}
    if (nargin >= 4) {blockSize  = args(3).int_value();}

    FloatNDArray input = args(0).float_array_value();

    int w = input.dim1();
    if (w != input.dim2()) {
      octave_stdout << "ERROR: Square matrix are required.";
      return octave_value_list();
    }
    dim_vector s = input.dims();
    
    FloatNDArray output(s);

    float *h_input  = (float*) input.fortran_vec();
    float *h_output = (float*)output.fortran_vec();

    haar(h_output, h_input, w, steps, inverse, blockSize);

    if (! error_state) {
      return octave_value (output);
    }
  }
  
  return octave_value_list();
}

