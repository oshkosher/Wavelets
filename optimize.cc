#include "optimize.h"


bool varyBinCounts(OptimizationData *optData,
                        float *thresholdValueOut, int *binCountOut) {

  int size;
  float l1, l2, mse, psnr, relErr;

  // Note: insert something smarter here

  /*
  // first try some different threshold settings without quantization
  printf("No quantization, test just threshold value.\n");
  for (float thresholdFrac = 0; thresholdFrac < .95; thresholdFrac += 0.1) {
    float thresh = 0;
    if (thresholdFrac > 0) {
      int offset = (int)(optData->count() * thresholdFrac);
      thresh = optData->getSorted(offset);
    }
    if (!testParameters(optData, thresh, -1, QUANT_ALG_UNKNOWN,
                        &size, &l1, &l2, &mse, &psnr)) return false;
    printf("  %2d%% (%.4f)  %.2f mse, %.2f psnr\n", (int)(thresholdFrac * 100),
           thresh, mse, psnr);
    fflush(stdout);
  }
  */

  float frac = .95;
  *thresholdValueOut = optData->getSorted((int)(optData->count() * frac));

  printf("\nthresh = %2.0f%%, log quant\n",
         frac * 100);
  for (int binCount=20; binCount <= 30000; binCount = (binCount*3)/2) {
    if (!testParameters(optData, *thresholdValueOut, binCount, QUANT_ALG_LOG,
                        &size, &l1, &l2, &mse, &psnr, &relErr)) return false;
    printf("  %2d bins, %d bytes, %.2f mse, %g relErr, %.2f psnr\n", binCount, size, mse, relErr, psnr);
    fflush(stdout);
  }    


  *thresholdValueOut = optData->getSorted(optData->count() / 2);
  *binCountOut = 100;

  return false;
}


bool fixedRelativeError(OptimizationData *optData, float *thresholdValueOut,
			int *binCountOut, float target) {
               

  float psnr, relativeError;
  int size;
  QuantizeAlgorithm quantAlg = QUANT_ALG_LOG;

  for (float frac = .15f; frac > .045; frac -= .01f) {
    float threshValue = optData->getSorted((int)(optData->count() * (1-frac)));
    printf("threshold %.1f%%", (1-frac) * 100);
    fflush(stdout);

    testParameters(optData, threshValue, 0, quantAlg, &size, 
                   NULL, NULL, NULL, &psnr, &relativeError);
    // printf(", max pSNR = %.2f, relErr = %.2e\n", psnr, relativeError);

    if (relativeError > target) {
      printf("  unable, minimum relErr = %.2e\n", relativeError);
      continue;
    }

    for (int binCount=50; binCount <= 30000; binCount = (int)(binCount*1.1)) {

      testParameters(optData, threshValue, binCount, quantAlg, &size, 
                     NULL, NULL, NULL, &psnr, &relativeError);
      
      if (relativeError < target) {
	printf("  %d bins: psnr=%.2f, relErr=%.2e, size=%d\n",
	       binCount, psnr, relativeError, size);
	fflush(stdout);
	break;
      }

    }

  }

  return false;
}

bool plotAll(OptimizationData *optData, float *thresholdValueOut,
             int *binCountOut) {

  float psnr;
  int size;
  QuantizeAlgorithm quantAlg = QUANT_ALG_LOG;

  printf("Threshold\tBins\tsize\tpSNR\n");

  for (float frac = .5f; frac < .999f; frac += 0.01f) {

    float threshValue = optData->getSorted((int)(optData->count() * frac));

    for (int binCount=100; binCount <= 2000; binCount = (int)(binCount*1.5)) {

      testParameters(optData, threshValue, binCount, quantAlg, &size, 
                     NULL, NULL, NULL, &psnr);
      
      printf("%.3f\t%d\t%d\t%.6g\n", frac, binCount, size, psnr);
      fflush(stdout);

    }
  }

  return false;
}


bool optimizeParameters(OptimizationData *optData,
                        float *thresholdValueOut, int *binCountOut) {
  // return varyBinCounts(optData, thresholdValueOut, binCountOut);

  // return fixedRelativeError(optData, thresholdValueOut, binCountOut, .0005);
  return plotAll(optData, thresholdValueOut, binCountOut);

}
