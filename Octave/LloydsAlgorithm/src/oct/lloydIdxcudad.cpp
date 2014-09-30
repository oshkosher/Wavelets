#include <octave/oct.h>

void setupLloyds(double *points, unsigned int psize, double pmax, double pmin);
void lloyd(double *codebook, unsigned int csize, double stop_criteria, double *table, double &dist, double &reldist, unsigned int *groups);
void finalized();

DEFUN_DLD (lloydIdxcudad, args, , "\n"
          "  lloydcudad(points, codebook, tol = 1e-7)\n"
          "  Do a lloyds quantization in double precision.\n"
          "\n") {
	octave_value_list retval;

  int nargin = args.length();

  if (nargin < 2) {
    print_usage();
  } else {

    NDArray points = args(0).array_value();

    unsigned int psize = points.dim1() * points.dim2();

    NDArray codebook = args(1).array_value();

    unsigned int csize = codebook.dim1() * codebook.dim2();

		double tol;

		if (nargin < 3) {
			tol = 1e-7;
		} else {
			tol = args(2).double_value();
		}

    double *h_points = (double*) points.fortran_vec();

		double maxpoints = h_points[0];
		double minpoints = h_points[0];
		for(int i = 1; i < psize; i++) {
			double aux = h_points[i];
			maxpoints = (aux > maxpoints) ? aux : maxpoints;
			minpoints = (aux < minpoints) ? aux : minpoints;;
		}

    uint32NDArray groups(points.dims());
		
		dim_vector dim(1);
		dim(0) = csize-1;
    NDArray table(dim);

    double *h_codebook       = (double*) codebook.fortran_vec();
    unsigned int *h_groups  = (unsigned int*) groups.fortran_vec();
    double *h_table          = (double*) table.fortran_vec();

		double dist;
		double reldist;

		double stop_criteria = abs(tol);

		setupLloyds(h_points, psize, maxpoints, minpoints);
		lloyd(h_codebook, csize, stop_criteria, h_table, dist, reldist, h_groups);
		finalized();

    if (! error_state) {
      retval(0) = octave_value(groups);
      retval(1) = octave_value(codebook);
    }
  }
  
  return retval;
}

