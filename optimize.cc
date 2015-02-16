#include "optimize.h"
#include "huffman.h"
#include "test_compress_common.h"
#include "test_compress_cpu.h"

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

    quantAlg - QUANT_ALG_LOG, QUANT_ALG_LLOYD, or QUANT_ALG_UNIFORM

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

  WaveletCompressionParam param;
  param.binCount = binCount;
  param.quantAlg = quantAlg;
  param.thresholdValue = thresholdValue;
  param.waveletAlg = o->waveletAlg;
  param.transformSteps = o->transformSteps;

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

    // quantize and dequantize the data

    CubeInt quantizedData;
    const float *nonzeroData = o->sorted + o->count - o->nonzeroCount;
    if (!quantize(*o->transformedData, quantizedData,
                  o->maxAbsVal, param,
                  nonzeroData, o->nonzeroCount,
                  o->minVal, o->maxVal))
      return false;
    
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


bool optimizeParameters(OptimizationData *optData,
                        float *thresholdValue, int *binCount) {

  int size;
  float l1, l2, mse, psnr;

  // Note: insert something smarter here

  // first try some different threshold settings without quantization
  printf("No quantization, test just threshold value.\n");
  for (float thresholdFrac = 0; thresholdFrac < .95; thresholdFrac += 0.1) {
    float thresh = 0;
    if (thresholdFrac > 0) {
      int offset = (int)(optData->count * thresholdFrac);
      thresh = optData->sorted[offset];
    }
    if (!testParameters(optData, thresh, -1, QUANT_ALG_UNKNOWN,
                        &size, &l1, &l2, &mse, &psnr)) return false;
    printf("  %2d%% (%.4f)  %.2f mse, %.2f psnr\n", (int)(thresholdFrac * 100),
           thresh, mse, psnr);
  }
  
  float frac = .7;
  *thresholdValue = optData->sorted[(int)(optData->count * frac)];

  printf("\nthresh = %2.0f%%, try different quantization bin counts\n",
         frac * 100);
  for (int binCount=10; binCount <= 2560; binCount *= 2) {
    if (!testParameters(optData, *thresholdValue, binCount, QUANT_ALG_LOG,
                        &size, &l1, &l2, &mse, &psnr)) return false;
    printf("  %2d bins, %d bytes, %.2f mse, %.2f psnr\n", binCount, size, mse, psnr);
  }    

  *binCount = 100;

  return false;
}
