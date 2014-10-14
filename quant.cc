#include "quant.h"


float quant_copysignf(float x, float s) {
  // for some reason, copysignf is defined on 64-bit Windows, but not 32-bit
#if defined(_WIN32) && !defined(_M_X64)
  if (s < 0)
    return -x;
  else
    return x;
#else
  return copysignf(x, s);
#endif
}


float quant_log2(float x) {
#ifdef _WIN32
  return (log(fabsf(x))/log(2.0));
#else
  return log2f(x);
#endif
}
