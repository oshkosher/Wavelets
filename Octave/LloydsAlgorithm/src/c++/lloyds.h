#define DEFAULT_LLOYD_STOP_CRITERIA 10e-7

void lloyd(const float *points, unsigned int psize,float *codebook, unsigned int csize, float *partition, float &dist, float &reldist, unsigned int *groups, float stop_criteria = DEFAULT_LLOYD_STOP_CRITERIA, int *iterations = NULL);
