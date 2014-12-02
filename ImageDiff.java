import java.awt.image.*;
import javax.imageio.*;
import java.io.*;

/*
  Read two data files and output a PNG image depicting where the
  data files differ in green and red pixels.
*/
public class ImageDiff {

  private static float findMaxDiff(WaveletSampleImage.Matrix data1,
                                   WaveletSampleImage.Matrix data2) {
    float max = 0;

    for (int i=0; i < data1.data.length; i++) {
      max = Math.max(max, Math.abs(data1.data[i] - data2.data[i]));
    }
    return max;
  }


  // Return the difference as an integer from 0 to 255, scaled
  // so that maxDiff returns 255.
  private static int scaleDiff(float diff, float maxDiff) {
    return (int)(diff / maxDiff * 255 + .5f);
  }
    

  public static void main(String[] args) throws Exception {
    
    if (args.length != 3) {
      System.out.println
        ("\n  ImageDiff <data1> <data2> <image_out>\n");
      return;
    }

    WaveletSampleImage.Matrix data1 = WaveletSampleImage.readMatrix(args[0]);
    WaveletSampleImage.Matrix data2 = WaveletSampleImage.readMatrix(args[1]);

    if (data1.width != data2.width ||
        data1.height != data2.height) {
      System.out.println("Data files are different sizes.");
      return;
    }

    float maxDiff = findMaxDiff(data1, data2);
    System.out.println("max diff = " + maxDiff);
    
    BufferedImage image = new BufferedImage
      (data1.width, data1.height, BufferedImage.TYPE_INT_RGB);

    int i = 0;
    for (int y=0; y < data1.height; y++) {
      for (int x=0; x < data1.width; x++) {
        float diff = data1.data[i] - data2.data[i];
        i++;

        int color = 0;

        if (diff > 0) {
          // red
          color = scaleDiff(diff, maxDiff) << 16;
        } else if (diff < 0) {
          // green
          color = scaleDiff(-diff, maxDiff) << 8;
        }

        image.setRGB(x, y, color);
      }
    }

    ImageIO.write(image, WaveletSampleImage.filenameSuffix(args[2]),
                  new File(args[2]));
  }
}

    
