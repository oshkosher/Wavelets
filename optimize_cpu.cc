#include "optimize.h"
#include "huffman.h"
#include "test_compress_common.h"
#include "test_compress_cpu.h"

float OptimizationData::getSorted(int index) {
  if (index < 0 || index >= count()) {
    assert(index >= 0 && index < count());
    return 0;
  }

  return sorted[index];
}


int OptimizationData::findSorted(float value) {
  return std::lower_bound(sorted, sorted+count(), value) - sorted;
}


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
bool testParameters(OptimizationData *o,
                    float thresholdValue, int binCount,
                    QuantizeAlgorithm quantAlg,
                    int *outputSizeBytes,
                    float *l1Error, float *l2Error, float *mse, float *pSNR) {

  WaveletCompressionParam param = o->transformedData->param;
  param.binCount = binCount;
  param.quantAlg = quantAlg;
  param.thresholdValue = thresholdValue;

  CubeFloat inverseWaveletInput;
  inverseWaveletInput.size = o->transformedData->size;
  inverseWaveletInput.allocate();

  *outputSizeBytes = 0;

  int count = o->transformedData->count();


  // If binCount is nonpositive, don't quantize and dequantize.
  // Just apply the threshold and reverse the wavelet transform.
  if (binCount <= 0) {

    // apply threshold
    const float *readp = o->transformedData->pointer(0,0,0);
    float *writep = inverseWaveletInput.pointer(0,0,0);

    for (int i=0; i < count; i++) {
      if (fabsf(readp[i]) <= thresholdValue) {
        writep[i] = 0;
      } else {
        writep[i] = readp[i];
      }
    }

  } else {  // binCount > 0

    // o->transformedData->print("before quantize");

    // quantize and dequantize the data
    CubeInt quantizedData;
    int firstNonzero = o->findSorted(thresholdValue);
    int nonzeroCount = count - firstNonzero;
    const float *nonzeroData = o->sorted + firstNonzero;
    if (!quantize(*o->transformedData, quantizedData,
                  o->maxAbsVal, param,
                  nonzeroData, nonzeroCount,
                  o->minVal, o->maxVal))
      return false;

    // quantizedData.print("after quantize");
    
    // compute the huffman coding
    Huffman huff;
    huff.init(binCount);
    int *readp = quantizedData.pointer(0,0,0);
    int *endp = readp + count;
    while (readp < endp)
      huff.increment(*readp++);
    huff.computeHuffmanCoding();

    // just get the size of the encoded data; don't write it out
    *outputSizeBytes = huff.totalEncodedLengthBytes();

    // dequantize
    if (!dequantize(quantizedData, inverseWaveletInput, param))
      return false;
  }

  // compare to original
  ErrorAccumulator errAccum;

  computeErrorRatesAfterDequant(&inverseWaveletInput, param, o->originalData,
                                &errAccum);

  *l1Error = errAccum.getL1Error();
  *l2Error = errAccum.getL2Error();
  *mse = errAccum.getMeanSquaredError();
  *pSNR = errAccum.getPeakSignalToNoiseRatio();

  return true;
}
