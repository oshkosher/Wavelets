#include <octave/oct.h>

void setupLloyds(float *points, unsigned int psize, unsigned int *groups, float pmax, float pmin);
void lloyd(float *codebook, unsigned int csize, float stop_criteria, float *table, float &dist, float &reldist);
void finalize();

DEFUN_DLD (lloydcuda, args, , "\n"
          "  lloydcuda(input, steps = 1, inverse = false, text = false, blockSize = -1)\n"
          "  Do a lloyds quantization.\n"
          "\n") {
	octave_value_list retval;

  int nargin = args.length();

  if (nargin < 2) {
    print_usage();
  } else {

    //if (nargin >= 2) {steps  = args(1).int_value();}
    //if (nargin >= 3) {inverse    = args(2).bool_value();}
    //if (nargin >= 4) {blockSize  = args(3).int_value();}

    FloatNDArray points = args(0).float_array_value();

    unsigned int psize = points.dim1() * points.dim2();

    FloatNDArray codebook = args(1).float_array_value();

    unsigned int csize = codebook.dim1() * codebook.dim2();

		float tol;

		if (nargin < 3) {
			tol = 1e-7;
		} else {
			tol = args(2).float_value();
		}

    float *h_points = (float*) points.fortran_vec();

		float maxpoints = h_points[0];
		float minpoints = h_points[0];
		for(int i = 1; i < psize; i++) {
			float aux = h_points[i];
			maxpoints = (aux > maxpoints) ? aux : maxpoints;
			minpoints = (aux < maxpoints) ? aux : maxpoints;;
		}

    uint32NDArray groups(psize);
		
    FloatNDArray table(csize);

    float *h_codebook       = (float*) codebook.fortran_vec();
    unsigned int *h_groups  = (unsigned int*) groups.fortran_vec();
    float *h_table          = (float*) table.fortran_vec();

		float dist;
		float reldist;

		float stop_criteria = abs(tol);

		setupLloyds(h_points, psize, h_groups, maxpoints, minpoints);
		lloyd(h_codebook, csize, stop_criteria, h_table, dist, reldist);
		finalize();

    if (! error_state) {
      retval(0) = octave_value(table);
      retval(1) = octave_value(codebook);
      retval(2) = octave_value(dist);
      retval(3) = octave_value(reldist);
    }

		retval(0) = octave_value(points);
  }
  
  return retval;
}

