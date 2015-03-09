#ifndef __OPTIMIZE_H__
#define __OPTIMIZE_H__

#include "wavelet.h"

class OptimizationData {
  friend bool testParameters
    (OptimizationData *o, float thresholdValue, int binCount,
     QuantizeAlgorithm quantAlg, int *outputSizeBytes,
     float *l1Error, float *l2Error, float *mse, float *pSNR,
     float *relativeError);
  
  const float *sorted;  // sorted absolute values, length is equal to count()

  // original cubelet data
  const Cube *originalData;

  // cubelet data after wavelet transform, 3-d array of floats
  const CubeFloat *transformedData;

  // const int nonzeroCount;  // # of values in sorted[] that are nonzero

  // compress strings of zeros in huffman encoding
  bool doCompressZeros;

 public:

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

  // if 3-d data, whether the wavelet transform was done as one 3-d tranform
  // or multiple 2-d transforms.
  const bool do2DTransform;
  
  OptimizationData(const Cube *originalData_,
                   const CubeFloat *transformedData_,
                   const float *sorted_, bool doCompressZeros_,
                   float minVal_, float maxVal_, float maxAbsVal_)
    :  sorted(sorted_), originalData(originalData_),
    transformedData(transformedData_),
    doCompressZeros(doCompressZeros_),
    minVal(minVal_), maxVal(maxVal_), maxAbsVal(maxAbsVal_),
    transformSteps(transformedData_->param.transformSteps),
    waveletAlg(transformedData_->param.waveletAlg),
    do2DTransform(transformedData_->param.do2DTransform) {}


  // returns the number of elements in the input data.
  int count() {return transformedData->count();}

  /**
     Return a value from the 'sorted' array, which has length 'count()'.
     This is a function because if the sorted array is on the GPU,
     the element will need to be copied from the GPU.
     Fails an assertion and return 0 if index is out of bounds.
  */
  float getSorted(int index);

  /**
    Return the index of the first value in 'sorted' that is greater than
    or equal to 'value'. Returns count() if 'value' is larger than every
    value in 'sorted'.
  */
  int findSorted(float value);
};
                    

/**
   Put the error control / parameter optimization logic here.
   
   This will be called by test_compress_cpu after the data is loaded
   and run through a wavelet transform. All the data is in "optData".
   You can read any of the public fields directly, but to access the
   'sorted' array ue the getSorted() and findSorted() accessor functions,
   because the sorted array might be on the GPU.

   This code in this function should call testParameters() to try out
   different settings for thresholdValue and binCount. When it is
   satisfied with the results, set the output parameters thresholdValue
   and binCount and return.

   Return true if you want the calling program to continue to process
   the data with these parameters and produce a compressed output file.
   
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

    outputSizeByte - size of the output data, in bytes
    l1Error, l2Error, pSNR, relativeSqError - error metrics

    On error, return false.
    On success, return true.
*/
bool testParameters(OptimizationData *optData,
                    float thresholdValue, int binCount,
                    QuantizeAlgorithm quantAlg,
                    int *outputSizeBytes,
                    float *l1Error = NULL, float *l2Error = NULL,
                    float *meanSquaredError = NULL,
                    float *pSNR = NULL, float *relativeError = NULL);

#endif // __OPTIMIZE_H__


