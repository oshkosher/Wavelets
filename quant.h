#ifndef __QUANT_H__
#define __QUANT_H__

//#include <cmath>
#include <math.h>
#include <vector>
#include <cassert>
#include <algorithm>
#include "nixtimer.h"

#ifdef __CUDACC__
#include "cuda.h"
#define HD __host__ __device__
#define DV __device__

#else

#define HD
#define DV

// For some reason, copysignf is defined on 64-bit Windows, but not 32-bit.
// And on 64-bit, it's called _copysignf, not copysignf.
#ifdef _WIN32
#ifdef _M_X64
#define copysignf _copysignf
#else
#define copysignf QuantFunction::quant_copysignf
#endif
#endif

#endif // __CUDACC__


class QuantFunction {
 public:
  HD static float quant_copysignf(float x, float s) {
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

  HD static float quant_log2(float x) {
#ifdef _WIN32
    return (log(fabsf(x))/log(2.0));
#else
    return log2f(x);
#endif
  }
};


/**
   Iterates through the data set calling the given quantizer object
   on each value.

   The quantizer object must implement two methods:
     int quant(float f);
     float dequant(int i);

*/
template <class Quantizer>
class QuantizationLooper {
  Quantizer *quantizer;
  double sumErrSquared, executeTime;
  size_t inputSize;

 public:
  QuantizationLooper(Quantizer *q = NULL) {
    init(q);
  }

  void init(Quantizer *q) {
    quantizer = q;
    executeTime = 0;
  }

  // XXX add peak signal-to-noise ratio
  void quantize(size_t length, const float *dataIn, int *dataOut = NULL,
		bool doComputeErr = false) {
    if (dataOut == NULL)
      dataOut = (int*) dataIn;

    sumErrSquared = 0;
    inputSize = length;

    double startTime = NixTimer::time();

    if (doComputeErr) {

      for (size_t i=0; i < length; i++) {
	float originalValue = dataIn[i];
	dataOut[i] = quantizer->quant(originalValue);
        // printf("%f\t%d\n", dataIn[i], dataOut[i]);
	float restoredValue = quantizer->dequant(dataOut[i]);
	// printf("%g\n", fabsf(originalValue - restoredValue));
	float err = restoredValue - originalValue;
	sumErrSquared += err*err;
      }

    } else {
      for (size_t i=0; i < length; i++) {
	dataOut[i] = quantizer->quant(dataIn[i]);
        // printf("%f\t%d\n", dataIn[i], dataOut[i]);
      }
    }

    executeTime = NixTimer::time() - startTime;
  }


  void dequantize(size_t length, int *dataIn, float *dataOut = NULL) {
    if (dataOut == NULL)
      dataOut = (float*) dataIn;

    double startTime = NixTimer::time();

    for (size_t i=0; i < length; i++) {
      dataOut[i] = quantizer->dequant(dataIn[i]);
    }

    executeTime = NixTimer::time() - startTime;
  }


  double getError() {
    return sumErrSquared / inputSize;
  }


  double getExecuteTime() {
    return executeTime;
  }

};


class QuantUniform {
  int bits, base;
  float threshold, maxVal;
  float scale, invScale;

 public:
  QuantUniform() {}

  QuantUniform(int bits_, float threshold_, float maxVal_) {
    init(bits_, threshold_, maxVal_);
  }

  HD void init(int bits_, float threshold_, float maxVal_) {
    bits = bits_;
    threshold = threshold_;
    maxVal = maxVal_;
    // if bits==8, then base = 127
    base = (1 << (bits-1)) - 1;

    /*
      y = base * ((x-threshold) / (maxVal-threshold))
        = (x-threshold) * (base / (maxVal-threshold))

      Let scale = base / (maxVal - threshold)
      and invScale = 1/scale

      y = (x-threshold) * scale
      y * invScale + threshold = x
    */
    
    scale = base / (maxVal - threshold);
    invScale = 1 / scale;

    printf("QuantUniform  bits=%d, threshold=%.8g, maxVal=%.8g, base=%d, scale = %.8g, invScale = %.8g\n", bits, threshold, maxVal, base, scale, invScale);
  }

  HD int quant(float x) const {

    float absx = fabsf(x);

    // apply threshold
    if (absx <= threshold) return 0;

    float scaled = (absx - threshold) * scale + 1;
    if (scaled >= base) scaled = base;

    return (int) copysignf(scaled, x);
  }

  HD float dequant(int x) const {
    if (x == 0) {
      return 0;
    } else if (x > 0) {
      return x * invScale + threshold;
    } else {
      return x * invScale - threshold;
    }
  }
};


class QuantLog {
  int bits, base;
  float threshold, invThresh, maxVal, lmax, lmaxInv, dqScale;

 public:
  QuantLog() {}
  
  QuantLog(int bits_, float threshold_, float maxVal_) {
    init(bits_, threshold_, maxVal_);
  }

  HD void init(int bits_, float threshold_, float maxVal_) {
    bits = bits_;
    threshold = threshold_;
    maxVal = maxVal_;

    base = (1 << (bits-1)) - 1;

    if (maxVal == threshold) {
      lmax = 1;
      lmaxInv = 1;
    } else {
      lmax = logf(maxVal/threshold);
      lmaxInv = 1 / lmax;
    }
    invThresh = (threshold == 0) ? 1 : (1 / threshold);
    dqScale = lmax / base;
  }

  HD int quant(float x) const {

    float absx = fabsf(x);
    if (absx <= threshold) return 0;
    // int sign=x/fabsf(x);
    
    float lnVal = logf(absx * invThresh);
    float result = base * lnVal * lmaxInv + 1;

    if (result > base) result = base;

    return (int) copysignf(result, x);
  }

  HD float dequant(int x) const {
    if (x == 0) return 0;
    // int sign=x/abs(x);
    float lnVal=fabsf(x*dqScale);
    return copysignf(threshold * expf(lnVal), x);
  }
};


class QuantCodebook {
  float lastBoundary;

 public:
  std::vector<float> boundaries, codebook;

  QuantCodebook() {}

  QuantCodebook(const std::vector<float> &boundaries_,
		const std::vector<float> &codebook_) {
    init(boundaries_, codebook_);
  }

  void init(const std::vector<float> &boundaries_,
	    const std::vector<float> &codebook_) {
    boundaries = boundaries_;
    codebook = codebook_;
    lastBoundary = boundaries[boundaries.size()-1];
  }


  // Generate boundaries and codebook entries based on bins with
  // equal numbers of values in each.
  void initCountBins(int count, float *data, int bits, float thresh);

  int quant(float x) const {

    // if the x is >= the end of the last bin, return the last bin
    if (x >= lastBoundary) return boundaries.size();
    /*
    int lastPos = (int)boundaries.size()-1;
    float lastBound = boundaries[lastPos];
    if (x >= lastBound)
      return lastPos+1;
    */

    // Find the first boundary that is greater than x.
    // For example, if the boundaries are:
    // boundaries[0] = 3
    // boundaries[1] = 5
    // boundaries[2] = 10
    // Given .5, it returns 0 because 3 > .5
    // Given 5, it returns 2, because 10 > 5
    // Given 100, it return 3, because no entry is > 100

    int bin = std::upper_bound(boundaries.begin(), boundaries.end(), x)
      - boundaries.begin();

    return bin;
  }

  float dequant(int x) const {
    assert(x >= 0 && x < (int)codebook.size());

    return codebook[x];
  }
  
};

#endif // __QUANT_H__
