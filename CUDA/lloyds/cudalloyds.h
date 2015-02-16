#define DEFAULT_LLOYD_STOP_CRITERIA 10e-7

void cudaLloyd(const float *points, unsigned int psize, float *codebook,
               unsigned int csize,
               float stop_criteria = DEFAULT_LLOYD_STOP_CRITERIA,
               bool points_are_on_gpu = false, int *iterations = NULL);
