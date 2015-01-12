import org.bytedeco.javacv.FFmpegFrameGrabber;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.awt.Graphics;
import java.io.File;

/*
  This requires JavaCV to compile or run. To download:
  http://search.maven.org/remotecontent?filepath=org/bytedeco/javacv/0.10/javacv-0.10-bin.zip

  Sample usage:
  java -cp '.;protobuf-2.6.0/protobuf.jar;c:/Apps/javacv-bin/javacv.jar' MovieToCubelets around.cube 2>/dev/null

  It produces excessive output on stderr, so you might
  want to add "2>/dev/null" to the command line.
*/

public class MovieToCubelets {

  public static final int CUBELET_WIDTH = 512;
  public static final int CUBELET_HEIGHT = 512;
  public static final int CUBELET_DEPTH = 512;

  public static void main(String[] args) throws Exception {
    if (args.length < 1 || args.length > 2) printHelp();

    String inputMovieFile = args[0];
    String outputCubeletFile = null;

    if (args.length > 1 && !args[1].equals("-"))
      outputCubeletFile = args[1];

    FFmpegFrameGrabber grabber = new FFmpegFrameGrabber(inputMovieFile);

    grabber.start();

    int width = grabber.getImageWidth();
    int height = grabber.getImageHeight();
    int frameCount = grabber.getLengthInFrames();

    if (args.length < 2) {
      System.out.printf("%dx%d, %.2f frames/sec, %d frames\n",
                        width, height, grabber.getFrameRate(), frameCount);
      grabber.stop();
      return;
    }

    CubeletFile.Writer cubeletOut = new CubeletFile.Writer();
    cubeletOut.open(outputCubeletFile);
    CubeletFile.Cubelet cube = new CubeletFile.Cubelet();
    cube.setSize(CUBELET_WIDTH, CUBELET_HEIGHT, CUBELET_DEPTH);
    cube.datatype = CubeletFile.Cubelet.DataType.UINT8;
    cube.byteData = new byte[CUBELET_WIDTH * CUBELET_HEIGHT * CUBELET_DEPTH];

    for (int frameNo=0; frameNo < CUBELET_DEPTH; frameNo++) {
      // get one frame
      BufferedImage image = grabber.grab().getBufferedImage();

      // convert to grayscale
      BufferedImage grayImage = new BufferedImage
        (width, height, BufferedImage.TYPE_BYTE_GRAY);
      Graphics g = grayImage.getGraphics();
      g.drawImage(image, 0, 0, null);
      g.dispose();
      image = grayImage;

      // crop one cubelet frame
      image = image.getSubimage(0, 0, CUBELET_WIDTH, CUBELET_HEIGHT);

      // extract the bytes
      getGrayscaleBytes(image, cube.byteData,
                        CUBELET_WIDTH * CUBELET_HEIGHT * frameNo);

      System.out.printf("\rFrame %d of %d", frameNo+1, CUBELET_DEPTH);
      System.out.flush();
    }
    System.out.println();

    grabber.stop();

    cubeletOut.addCubelet(cube);
    cubeletOut.close();
  }

  
  public static void getGrayscaleBytes
    (BufferedImage image, byte[] byteData, int offset) {
    
    int width = image.getWidth(), height = image.getHeight();
    for (int y=0; y < height; y++) {
      for (int x=0; x < width; x++) {
        byteData[offset++] = (byte) (image.getRGB(x, y) & 0xff);
      }
    }
  }
  

  static void printHelp() {
    System.out.println("\n  MovieToCubelets <input_movie>\n" +
                       "    Output basic settings about the given movie.\n" +
                       "  MovieToCubelets <input_movie> <output_cubelet_file>\n" +
                       "    Convert sequences of frames into cubelets.\n" +
                       "    Set output_cubelet_file to '-' for stdout.\n");
    System.exit(1);
  }
}
