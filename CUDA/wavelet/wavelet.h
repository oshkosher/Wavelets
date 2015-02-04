#define KERNEL_RADIUS 4
#define KERNEL_LENGTH (2 * KERNEL_RADIUS + 1)

extern "C" void setUpFilter(const float *filter);

extern "C" void fwt_1D_GPU(float *data, const unsigned level, const unsigned nx, const unsigned ny);
extern "C" void iwt_1D_GPU(float *data, const unsigned level, const unsigned nx, const unsigned ny);

extern "C" void  comp(float *data, const unsigned nx, const unsigned ny, const unsigned nz, const unsigned lvlx, const unsigned lvly, const unsigned lvlz);
extern "C" void ucomp(float *data, const unsigned nx, const unsigned ny, const unsigned nz, const unsigned lvlx, const unsigned lvly, const unsigned lvlz);
