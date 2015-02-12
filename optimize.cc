#include "optimize.h"


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
