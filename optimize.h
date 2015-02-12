#ifndef __OPTIMIZE_H__
#define __OPTIMIZE_H__

#include "wavelet.h"

struct OptimizationData {

  // original cubelet data, 3-d array of unsigned char
  const CubeByte *originalData;

  // cubelet data after wavelet transform, 3-d array of floats
  const CubeFloat *transformedData;

  const int count;  // # of values in sorted[]

  const float *sorted;  // sorted absolute values

  const int nonzeroCount;  // # of values in sorted[] that are nonzero

  const float minVal;  // minimum value in the data (will be negative)

  const float maxVal;  // maximum value in the data

  const float maxAbsVal;  // maximum absolute value in the data, should be
                          // max(abs(minVal), abs(maxVal))

  // # of transform steps along each axis.
  // This is a structure of three integers, defined in wavelet.h:
  // transformSteps.x, transformSteps.y, transformSteps.z
  const scu_wavelet::int3 transformSteps;

  // Wavelet algorithm used, WAVELET_CDF97 by default.
  const WaveletAlgorithm waveletAlg;


  OptimizationData(const CubeByte *od, const CubeFloat *td,
                   int n, const float *sd, int nzn,
                   float min, float max, float amax,
                   scu_wavelet::int3 xfs, WaveletAlgorithm wa)
  : originalData(od), transformedData(td), count(n), sorted(sd),
    nonzeroCount(nzn), 
    minVal(min), maxVal(max), maxAbsVal(amax),
    transformSteps(xfs), waveletAlg(wa) {}
};
                    

/**
   Put the error control / parameter optimization logic here.
   
   This will be called by test_compress_cpu after the data is loaded
   and run through a wavelet transform. All the data is in "optData".
   See the definition of "struct OptimizationData" above. Consider the
   data read-only, but there is data in it that may be useful, particularly
   the "sorted" array.

   This code in this function should call testParameters() to try out
   different settings for thresholdValue and binCount. When it is
   satisfied with the results, set the output parameters thresholdValue
   and binCount and return.

   Return true if you want the calling program to continue to process
   the data with these parameters and produced a compressed output file.
   
   Return false if this was just an experiment and you want the calling
   program to exit without producing a compressed output file.
 */
bool optimizeParameters
(OptimizationData *optData, // input data, after being wavelet-transformed

                          // output parameters:
 float *thresholdValue,   // threshold value
 int *binCount);          // quantization bin count



/** Call this to try out compression parameters thresholdValue and binCount.
    The resulting output file size and some error metrics will be returned.

    input parameters:

    optData - just pass the value passed to optimizeParameters()
    
    thresholdValue - the threshold value to test (the actual value, not
                     the fraction of the data set that is to be cut off)

    binCount - the quantization bin count value to test.
               If this is negative, no quantization will be done.
               testParameters() will just apply the threshold and
               then the inverse wavelet transform.

    quantAlg - QUANT_ALG_LOG, QUANT_ALG_LLOYD, QUANT_ALG_UNIFORM,
               or QUANT_ALG_UNKNOWN if no quantization is to be done.

    output parameters:

    outputSize - size of the output data, in bytes
    l1Error, l2Error, pSNR - error metrics

    On error, return false.
    On success, return true.
*/
bool testParameters(OptimizationData *optData,
                    float thresholdValue, int binCount,
                    QuantizeAlgorithm quantAlg,
                    int *outputSizeBytes,
                    float *l1Error, float *l2Error, float *meanSquaredError,
                    float *pSNR);

#endif // __OPTIMIZE_H__


