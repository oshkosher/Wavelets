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

  For FFmpegFrameGrabber, see the source:
    https://code.google.com/p/javacv/source/browse/src/main/java/com/googlecode/javacv/FFmpegFrameGrabber.java?name=
*/

public class MovieToCubelets {

  public static final int CUBELET_WIDTH = 512;
  public static final int CUBELET_HEIGHT = 512;
  public static final int CUBELET_DEPTH = 512;

  
  private static int cubesAcross, cubesDown;

  public static void main(String[] args) throws Exception {
    if (args.length < 1 || args.length > 2) printHelp();

    String inputMovieFile = args[0];
    String outputCubeletFile = null;

    if (args.length > 1 && !args[1].equals("-"))
      outputCubeletFile = args[1];

    // start up the frame grabber
    FFmpegFrameGrabber grabber = new FFmpegFrameGrabber(inputMovieFile);
    grabber.start();

    // get the size of the movie
    int width = grabber.getImageWidth();
    int height = grabber.getImageHeight();
    int frameCount = grabber.getLengthInFrames();

    // if an output file was not specified, just print the movie size
    if (args.length < 2) {
      System.out.printf("%dx%d, %.2f frames/sec, %d frames\n",
                        width, height, grabber.getFrameRate(), frameCount);
      grabber.stop();
      return;
    }

    // start up the cubelet output file
    CubeletFile.Writer cubeletOut = new CubeletFile.Writer();
    cubeletOut.width = width;
    cubeletOut.height = height;
    cubeletOut.open(outputCubeletFile);
    CubeletFile.Cubelet cubes[] = setupCubelets(width, height);

    int frameNo, currentDepth = 0, endFrame = CUBELET_DEPTH*3/2;
    for (frameNo=0; frameNo < endFrame; frameNo++) {
      // get one frame
      BufferedImage image = grabber.grab().getBufferedImage();

      // convert to grayscale
      BufferedImage grayImage = new BufferedImage
        (width, height, BufferedImage.TYPE_BYTE_GRAY);
      Graphics g = grayImage.getGraphics();
      g.drawImage(image, 0, 0, null);
      g.dispose();
      image = grayImage;

      // extract rectangles
      int cubeNo = 0;
      for (int y = 0; y < cubesDown; y++) {
        for (int x = 0; x < cubesAcross; x++) {
          CubeletFile.Cubelet cube = cubes[cubeNo++];
          BufferedImage subImage = image.getSubimage
            (cube.xOffset, cube.yOffset, cube.width, cube.height);
          getGrayscaleBytes(subImage, cube.byteData,
                            cube.width * cube.height * currentDepth);
          cube.zOffset = frameNo;
        }
      }

      // when the cubelets fill up, write them out and reset
      currentDepth++;
      if (currentDepth == CUBELET_DEPTH) {
        System.out.println();
        for (CubeletFile.Cubelet cube : cubes) {
          cube.depth = currentDepth;
          cubeletOut.addCubelet(cube);
        }
        currentDepth = 0;
      }

      System.out.printf("\rFrame %d of %d", frameNo+1, endFrame);
      System.out.flush();
    }
    System.out.println();

    grabber.stop();

    cubeletOut.depth = frameNo;

    // write out leftover frames
    if (currentDepth > 0) {
      for (CubeletFile.Cubelet cube : cubes) {
        cube.depth = currentDepth;
        cubeletOut.addCubelet(cube);
      }
    }

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


  /**
     Set up multple simultaneous cubelets to which data will be output.
     For example, if the movie is 150x100 and our cubelets are 64x64x64,
     then each frame will have 6 rectangles of data:
        across: 0..63, 64..127, 128..149
        down:   0..63, 64..99
  */
  static CubeletFile.Cubelet[] setupCubelets(int width, int height) {
  
    cubesAcross = (width + CUBELET_WIDTH - 1) / CUBELET_WIDTH;
    cubesDown = (height + CUBELET_HEIGHT - 1) / CUBELET_HEIGHT;

    CubeletFile.Cubelet cubes[] = 
      new CubeletFile.Cubelet[cubesAcross * cubesDown];
    
    int cubeNo = 0;
    for (int y = 0; y < cubesDown; y++) {
      for (int x = 0; x < cubesAcross; x++) {
        CubeletFile.Cubelet cube = new CubeletFile.Cubelet();
        cube.setSize(CUBELET_WIDTH, CUBELET_HEIGHT, CUBELET_DEPTH);
        cube.setOffset(x * CUBELET_WIDTH, y * CUBELET_HEIGHT, 0);
        cube.datatype = CubeletFile.Cubelet.DataType.UINT8;

        if (cube.yOffset + CUBELET_HEIGHT > height)
          cube.height = height - cube.yOffset;

        if (cube.xOffset + CUBELET_WIDTH > width)
          cube.width = width - cube.xOffset;

        cube.byteData = new byte[cube.width * cube.height * cube.depth];

        cubes[cubeNo++] = cube;
      }
    }

    return cubes;
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
