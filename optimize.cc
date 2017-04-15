#include "optimize.h"


bool varyThreshold(OptimizationData *optData,
                   float *thresholdValueOut, int *binCountOut) {
  int size;
  float l1, l2, mse, psnr, relErr;

  for (float invThresh = .5f; invThresh > .01; invThresh *= .8) {
    float thresh = 0;
    int offset = (int)(optData->count() * (1-invThresh));
    thresh = optData->getSorted(offset);

    if (!testParameters(optData, thresh, -1, QUANT_ALG_UNKNOWN,
                        &size, &l1, &l2, &mse, &psnr)) return false;
    printf("  %6.3f%% (%7.4f)  %5.2f mse, %5.2f psnr\n", (1-invThresh) * 100,
           thresh, mse, psnr);
    fflush(stdout);
    *thresholdValueOut = thresh;
  }

  *binCountOut = 256;
  return true;
}


bool varyBinCounts(OptimizationData *optData,
                   float *thresholdValueOut, int *binCountOut) {

  int size;
  float l1, l2, mse, psnr, relErr;

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


bool fixedPSNR(OptimizationData *optData, float *thresholdValueOut,
               int *binCountOut, float targetPSNR) {
               

  float psnr;
  int size;
  QuantizeAlgorithm quantAlg = QUANT_ALG_LLOYD;

  for (float frac = .40f; frac > .04f; frac -= .05) {
    float threshValue = optData->getSorted((int)(optData->count() * (1-frac)));
    printf("threshold %.3f", 1-frac);

    testParameters(optData, threshValue, 0, quantAlg, &size, 
                   NULL, NULL, NULL, &psnr);
    printf(", max pSNR = %.2f\n", psnr);
    fflush(stdout);
    if (psnr < targetPSNR) {
      continue;
    }

    bool found = false;
    for (int binCount=50; binCount <= 30000; binCount = (int)(binCount*1.1)) {

      testParameters(optData, threshValue, binCount, quantAlg, &size, 
                     NULL, NULL, NULL, &psnr);
      
        printf("  %d bins: psnr=%.2f, size=%d\n", binCount, psnr, size);
        fflush(stdout);
      if (psnr > targetPSNR) {
        found = true;
        break;
      }

    }
    if (!found) printf("\n");
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

    for (int binCount=100; binCount <= 2000; binCount = (int)(binCount*1.1)) {

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
  // return varyThreshold(optData, thresholdValueOut, binCountOut);

  // return varyBinCounts(optData, thresholdValueOut, binCountOut);

  // return fixedPSNR(optData, thresholdValueOut, binCountOut, 46);

  return plotAll(optData, thresholdValueOut, binCountOut);

}
