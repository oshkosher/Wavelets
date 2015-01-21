#define KERNEL_RADIUS 4
#define KERNEL_LENGTH (2 * KERNEL_RADIUS + 1)

extern "C" void setUpFilter(const float *filter);

extern "C" void    comp(float *data, const unsigned nx, const unsigned ny, const unsigned nz, const unsigned lvlx, const unsigned lvly, const unsigned lvlz);
extern "C" void invComp(float *data, const unsigned nx, const unsigned ny, const unsigned nz, const unsigned lvlx, const unsigned lvly, const unsigned lvlz);
