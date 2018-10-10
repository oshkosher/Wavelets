import org.bytedeco.javacv.FFmpegFrameGrabber;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;

/*
  This requires javacv.jar to compile or run, and 
     c:/Apps/javacv-bin/javacv.jar
  http://search.maven.org/remotecontent?filepath=org/bytedeco/javacv/0.10/javacv-0.10-bin.zip

  javac -cp c:/Apps/javacv-bin/javacv.jar GrabFrame.java
  java -cp '.;c:/Apps/javacv-bin/javacv.jar' GrabFrame LargeImages/AroundTheLake.MP4
*/

public class GrabFrame {

  public static void main(String[] args) throws Exception {
    if (args.length != 1) {
      System.out.println("\n  GrabFrame <filename>\n");
      return;
    }

    FFmpegFrameGrabber g = new FFmpegFrameGrabber(args[0]);

    g.start();

    System.out.printf("%dx%d, %.2f frames/sec, %d frames\n",
                      g.getImageWidth(), g.getImageHeight(),
                      g.getFrameRate(), g.getLengthInFrames());

    /*
    for (int i=0; i < 240; i++) {
      BufferedImage img = g.grab().getBufferedImage();
    }
    */

    for (int i=0; i < 1; i++) {
      BufferedImage img = g.grab().getBufferedImage();
      File filename = new File("frame"+i+".png");
      ImageIO.write(img, "png", filename);
      System.out.println("Wrote " + filename);
    }

    g.stop();
  }
}
