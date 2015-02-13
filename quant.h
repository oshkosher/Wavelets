#ifndef __QUANT_H__
#define __QUANT_H__

//#include <cmath>
#include <math.h>
#include <vector>
#include <cassert>
#include <algorithm>
#include "nixtimer.h"
#include "cucheck.h"

#ifndef __CUDACC__

// For some reason, copysignf is defined on 64-bit Windows, but not 32-bit.
// And on 64-bit, it's called _copysignf, not copysignf.
#ifdef _WIN32
#ifdef _M_X64
#define copysignf _copysignf
#else
#define copysignf Quantizer::quant_copysignf
#endif
#endif

#endif // __CUDACC__


class Quantizer {
  
 public:
  virtual void quantizeRow(const float *in, int *out, int count) = 0;
  virtual void dequantizeRow(const int *in, float *out, int count) = 0;
  int quant(float f) {
    int result;
    quantizeRow(&f, &result, 1);
    return result;
  }

  virtual ~Quantizer() {}

  HD static float copysignf(float x, float s) {
    // for some reason, copysignf is defined on 64-bit Windows, but not 32-bit
#if defined(_WIN32) && !defined(_M_X64)
    if (s < 0)
      return -x;
    else
      return x;
#else
    return ::copysignf(x, s);
#endif
  }

  HD static float log2(float x) {
#ifdef _WIN32
    return (logf(fabsf(x))/logf(2.0));
#else
    return ::log2f(x);
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
  int binCount;

 public:
  QuantizationLooper() : quantizer(NULL), binCount(0) {}

  QuantizationLooper(Quantizer *q, int b) {
    init(q, b);
  }

  void init(Quantizer *q, int b) {
    quantizer = q;
    binCount = b;
    executeTime = 0;
    sumErrSquared = 0;
  }

  // XXX add peak signal-to-noise ratio
  void quantize(size_t length, const float *dataIn, int *dataOut = NULL,
		bool doComputeErr = false) {
    assert(quantizer != NULL && binCount > 1);

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
        assert(dataOut[i] >= 0 && dataOut[i] < binCount);
      }

    } else {
      for (size_t i=0; i < length; i++) {
	dataOut[i] = quantizer->quant(dataIn[i]);
        // printf("%f\t%d\n", dataIn[i], dataOut[i]);
        assert(dataOut[i] >= 0 && dataOut[i] < binCount);
      }
    }

    executeTime = NixTimer::time() - startTime;
  }


  void dequantize(size_t length, const int *dataIn, float *dataOut = NULL) {
    assert(quantizer != NULL && binCount > 1);

    if (dataOut == NULL)
      dataOut = (float*) dataIn;

    double startTime = NixTimer::time();

    for (size_t i=0; i < length; i++) {
      assert(dataIn[i] >= 0 && dataIn[i] < binCount);
      dataOut[i] = quantizer->dequant(dataIn[i]);
      // printf("%d\t%f\n", dataIn[i], dataOut[i]);
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


/**
   Quantize uniformly over a specified range.

   Values < -maxVal map to 0
   Values > +maxVal map to binCount-1

   zeroBin = floor( (binCount-1) / 2 )

   For example, if there are 7 bins, zeroBin = 3:
     0..2 : negative values,
     3 : -threshold .. +threshold
     4..6 : positive values

   If there is an odd number of bins, then the number of bins for
   negative and positive values are the same. If even, then we'll give one
   more to the positives.

   For example, if there are 8 bins:
     0..2 : negative values,
     3 : -threshold .. +threshold
     4..7 : positive values
*/
class QuantUniform {
  int binCount;
  float threshold, maxVal, negOffset, posOffset;
  float negScale, negInvScale, posScale, posInvScale;

 public:

  HD void init(int binCount_, float threshold_, float maxVal_) {
    binCount = binCount_;
    threshold = threshold_;
    maxVal = maxVal_;

    assert(maxVal > threshold);
    assert(threshold >= 0);
    assert(binCount > 0);

    /*
      For each side (negative and positive)
      Map thresh..max to 1..binCountThisSide
        y = x * scale + offset
        1 = thresh * scale + offset
        N = max * scale + offset

        offset = 1 - thresh * scale
        offset = N - max * scale

        1 - thresh * scale = N - max * scale
        (max - thresh) * scale = N - 1

        scale = (N - 1) / (max - thresh)
        offset = 1 - thresh * scale

        y = x * scale + offset
        y - offset = x * scale
        (y - offset) * (1/scale) = x
        cache 1/scale as invScale
    */

    // number of bins for negative,positive numbers
    int negBinCount = (binCount-1) / 2;
    negScale = (negBinCount - 1) / (maxVal - threshold);
    

    int posBinCount = binCount / 2;
    posScale = (posBinCount - 1) / (maxVal - threshold);
    
    negInvScale = negScale==0 ? 0 : 1 / negScale;
    posInvScale = posScale==0 ? 0 : 1 / posScale;

    negOffset = 1 - threshold * negScale;
    posOffset = 1 - threshold * posScale;

    // printf("QuantUniform  bits=%d, threshold=%.8g, maxVal=%.8g, base=%d, scale = %.8g, invScale = %.8g\n", bits, threshold, maxVal, base, scale, invScale);
  }

  HD int quant(float x) const {

    int zeroBin = (binCount-1) >> 1;

    if (x < 0) {

      if (x >= -threshold) {
        return zeroBin;
      } else {
        x = x * negScale - negOffset;
        int bin = (int)x + zeroBin;
        if (bin < 0) {
          return 0;
        } else {
          return bin;
        }
      }

    } else {  // x >= 0

      if (x <= threshold) {
        return zeroBin;
      } else {
        x = x * posScale + posOffset + zeroBin;
        if (x > binCount - 1) {
          return binCount - 1;
        } else {
          return (int) x;
        }
      }
    }

  }

  HD float dequant(int x) const {
    assert(x >= 0 && x < binCount);

    int zeroBin = (binCount-1) >> 1;

    // add half a bin to dequantize into the midpoint of the range
    // of values that quantized into that bin
    if (x < zeroBin) 
      return (x - 0.5f + negOffset - zeroBin) * negInvScale;
    else if (x > zeroBin)
      return (x + 0.5f - posOffset - zeroBin) * posInvScale;
    else
      return 0;
  }

};


// Provide the services of a QuantUniform but with a virtual methods
class QuantizerUniform : public Quantizer {
  QuantUniform q;
  
 public:
  QuantizerUniform(int binCount, float threshold, float maxVal) {
    q.init(binCount, threshold, maxVal);
  }

  // Override
  virtual void quantizeRow(const float *in, int *out, int count) {
    for (int i=0; i < count; i++) out[i] = q.quant(in[i]);
  }

  // Override
  virtual void dequantizeRow(const int *in, float *out, int count) {
    for (int i=0; i < count; i++) out[i] = q.dequant(in[i]);
  }
};


/**
   Similar to QuantUniform:

   zeroBin = floor( (binCount-1) / 2 )

     0..(zeroBin-1) : negative
     zeroBin : -threshold .. +threshold
     (zeroBin+1)..(binCount-1) : positive
*/
class QuantLog {
  int binCount;
  float threshold;
  float negScale, negInvScale, negOffset;
  float posScale, posInvScale, posOffset;

 public:
  // this object is used as a CUDA shared memory variable, so it cannot
  // have any constructors
  HD QuantLog() {}

  HD void init(int binCount_, float threshold_, float maxVal_) {
    binCount = binCount_;
    threshold = threshold_;
    float maxVal = maxVal_;

    assert(maxVal > threshold);
    assert(threshold > 0);
    assert(binCount > 0);

    /*
      Map lo..hi to 1..N logarithmically
      y = 1 + scale * log(x / thresh)
        = 1 + scale * (log(x) - log(thresh))
        = 1 - scale*log(thresh) + scale * log(x)

      N = 1 + scale * log(maxVal / thresh)
      scale = (N - 1) / log(maxVal / thresh)
      offset = 1 - scale*log(thresh)
      y = scale * log(x) + offset

      (y - offset) / scale = log(x)

      e^( (y - offset) / scale ) = x
      x = e^( (y - offset) * invScale )
    */


    int negBinCount = (binCount-1) / 2;
    negScale = (negBinCount - 1) / logf(maxVal / threshold);
    negInvScale = 1 / negScale;
    negOffset = 1 - negScale * logf(threshold);

    int posBinCount = binCount / 2;
    posScale = (posBinCount - 1) / logf(maxVal / threshold);
    posInvScale = 1 / posScale;
    posOffset = 1 - posScale * logf(threshold);
  }

  HD int quant(float x) const {

    int zeroBin = (binCount-1) >> 1;

    if (x < 0) {

      x = -x;

      if (x <= threshold) {
        return zeroBin;
      } else {
        x = zeroBin - negScale * logf(x) - negOffset;
        if (x < 0) {
          return 0;
        } else {
          return (int) x;
        }
      }

    } else {  // x >= 0

      if (x <= threshold) {
        return zeroBin;
      } else {
        x = zeroBin + posScale * logf(x) + posOffset;
        if (x >= binCount) {
          return binCount-1;
        } else {
          return (int) x;
        }
      }
    }

  }

  HD float dequant(int x) const {
    assert(x >= 0 && x < binCount);

    int zeroBin = (binCount-1) >> 1;
    float sign, scale, offset;
    if (x < zeroBin) {
      sign = -1;
      offset = negOffset;
      scale = negInvScale;
    } else {
      sign = 1;
      offset = posOffset;
      scale = posInvScale;
    }

    if (x == zeroBin) {
      return 0;
    } else {
      return sign * expf( (sign*(x - zeroBin) - offset) * scale );
    }

      /*
    if (x < zeroBin)
      return -expf( (zeroBin - x - negOffset) * negInvScale );
    else if (x > zeroBin)
      return expf( (x - zeroBin - posOffset) * posInvScale );
    else
      return 0;
      */
  }
};


// Provide the services of a QuantLog but with a virtual methods
class QuantizerLog : public Quantizer {
  QuantLog q;
  
 public:
  QuantizerLog(int binCount, float threshold, float maxVal) {
    q.init(binCount, threshold, maxVal);
  }

  // Override
  virtual void quantizeRow(const float *in, int *out, int count) {
    for (int i=0; i < count; i++) out[i] = q.quant(in[i]);
  }

  // Override
  virtual void dequantizeRow(const int *in, float *out, int count) {
    for (int i=0; i < count; i++) {
      int tmp = in[i];
      out[i] = q.dequant(tmp);
      // printf("dequant %6d: %d %.2g\n", i, tmp, out[i]);
    }
  }
};


class QuantCodebook {
  float lastBoundary;

 public:
  std::vector<float> boundaries, codebook;

  /*
  QuantCodebook() {}

  QuantCodebook(const std::vector<float> &codebook_) {
    init(codebook_);
  }

  QuantCodebook(const std::vector<float> &boundaries_,
		const std::vector<float> &codebook_) {
    init(boundaries_, codebook_);
  }
  */

  // given both a codebook and set of boundaries
  void init(const std::vector<float> &codebook_,
            const std::vector<float> &boundaries_) {
    codebook = codebook_;
    boundaries = boundaries_;
    lastBoundary = boundaries[boundaries.size()-1];
  }

  // Given just the codebook, automatically create boundaries at the
  // midpoint between each pair of adjacent codebook entries.
  /*
  void init(const std::vector<float> &codebook_) {
    codebook = codebook_;
    boundaries.resize(codebook.size() - 1);

    for (size_t i = 0; i < boundaries.size(); i++) {
      boundaries[i] = (codebook[i] + codebook[i+1]) / 2;
    }

    lastBoundary = boundaries[boundaries.size()-1];
  }
  */

  // Generate boundaries and codebook entries based on bins with
  // equal numbers of values in each.
  void initCountBins(int count, float *data, int bits, float thresh);

  void printCodebook();

  int quant(float x) const {

    // if the x is >= the end of the last bin, return the last bin
    if (x >= lastBoundary) return (int)boundaries.size();

    int bin = (int)(std::lower_bound(boundaries.begin(), boundaries.end(), x)
      - boundaries.begin());

    return bin;

    // return quant(x, (int)boundaries.size(), boundaries.data());
  }

  static HD int quant(float x, int boundaryCount, const float *boundaries) {

    int count, step;
    const float *start = boundaries, *it;
    count = boundaryCount;

    // binary search, equivalent to std::lower_bound
    while (count > 0) {
      it = start;
      step = count / 2;
      it += step;
      if (*it < x) {
        start = it + 1;
        count -= step + 1;
      } else {
        count = step;
      }
    }

    return start - boundaries;
  }


  float dequant(int x) const {
    assert(x >= 0 && x < (int)codebook.size());

    return codebook[x];
  }

};


// Provide the services of a QuantCodebook but with virtual methods
class QuantizerCodebook : public Quantizer {
  QuantCodebook q;
  
 public:
  QuantizerCodebook(const std::vector<float> &codebook,
                    const std::vector<float> &boundaries) {
    q.init(codebook, boundaries);
  }

  // Override
  virtual void quantizeRow(const float *in, int *out, int count) {
    for (int i=0; i < count; i++) out[i] = q.quant(in[i]);
  }

  // Override
  virtual void dequantizeRow(const int *in, float *out, int count) {
    for (int i=0; i < count; i++) out[i] = q.dequant(in[i]);
  }
};

#endif // __QUANT_H__
