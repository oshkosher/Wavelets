/*
  Convert between JPEG images and 2-d arrays of grayscale values.

  Given a color image of any size:
    1. Crop the long direction evenly on both ends to make the image square.
    2. Resize to the next smaller power of 2.
    3. Convert to grayscale.
    4. Output as 2-d array of floats from 0 to 1.

  Given a 2-d array of floats from 0 to 1:
    1. Convert each value to a grayscale pixel.
    2. Output as a JPEG.

  Binary data format is described in data_io.h

  Ed Karrels, June 2014
*/

import java.awt.*;
import java.awt.color.*;
import java.awt.geom.*;
import java.awt.image.*;
import java.io.*;
import java.nio.file.*;
import java.text.*;
import java.util.*;
import javax.imageio.*;


public class WaveletSampleImage {

  private static final boolean VERBOSE_OUTPUT = false;

  static class Matrix {
    float[] data;
    int width, height;
  }


  public static void main(String args[]) {
    if (args.length < 1) printHelp();

    boolean textOutput = false;
    int argNo = 0;
    while (argNo < args.length && args[argNo].charAt(0)=='-') {
      if (args[argNo].equals("-text")) {
        textOutput = true;
        argNo++;
      }

      else printHelp();
    }

    // no input file listed
    if (argNo >= args.length) printHelp();

    String inputFile = args[argNo++], outputFile = null;
    if (argNo < args.length) outputFile = args[argNo++];


    String suffix = filenameSuffix(inputFile).toLowerCase();
    boolean sourceIsImage;

    if (suffix.equals("jpg") || suffix.equals("jpeg")) {
      sourceIsImage = true;
    } else if (suffix.equals("data")) {
      sourceIsImage = false;
    } else {
      System.err.println("Source file unrecognized. Should be .jpg, .jpeg, or .data.");
      return;
    }

    if (outputFile != null) {
      String destExt = filenameSuffix(outputFile).toLowerCase();
      if (sourceIsImage) {
        if (!destExt.equals("data")) {
          System.err.println("Output file must have .data extension.");
          return;
        }
      } else {
        if (!destExt.equals("jpg") && !destExt.equals("jpeg")) {
          System.err.println("Output file must have .jpg or .jpeg extension.");
          return;
        }
      }
    } else {
      String destExt = sourceIsImage ? "data" : "jpg";
      outputFile = inputFile.substring
        (0, inputFile.length() - suffix.length()) + destExt;
    }

    if (!new File(inputFile).exists()) {
      System.err.println("\"" + inputFile + "\" not found.");
      return;
    }

    if (sourceIsImage) {
      imageToMatrix(inputFile, outputFile, textOutput);
    } else {
      matrixToImage(inputFile, outputFile);
    }
  }

  public static void printHelp() {
    System.out.println
      ("\n" +
       "  WaveletSampleImage [-text] <src> [<dest>]\n" +
       "  If src is a jpeg image, convert it to grayscale data, cropped and\n" +
       "    resized so the width and height are equal and a power of two.\n" +
       "  If src is data, convert the float data into a grayscale image\n" +
       "  -text : output in text format rather than binary.\n"
       );
    System.exit(0);
  }

  public static String filenameSuffix(String str) {
    int dot = str.lastIndexOf(".");
    if (dot == -1) return "";
    return str.substring(dot+1);
  }


  /* Reads a JPEG image, convert to grayscale text data. */
  public static void imageToMatrix(String inputFile, String outputFile,
                                   boolean textOutput) {
    if (VERBOSE_OUTPUT) {
      System.out.println("Read image");
      System.out.flush();
    }
    BufferedImage image;
    try {
      image = ImageIO.read(new File(inputFile));
    } catch (IOException e) {
      System.err.println("Error reading \"" + inputFile + "\": " + e);
      return;
    }

    // crop image to be square
    if (VERBOSE_OUTPUT) {
      System.out.println("Crop");
      System.out.flush();
    }
    int size = image.getWidth();
    if (image.getWidth() > image.getHeight()) {
      size = image.getHeight();
      int leftCrop = (image.getWidth() - image.getHeight()) / 2;

      image = image.getSubimage(leftCrop, 0, image.getHeight(),
                                image.getHeight());

    } else if (image.getHeight() > image.getWidth()) {
      int topCrop = (image.getHeight() - image.getWidth()) / 2;
      image = image.getSubimage(0, topCrop, image.getWidth(),
                                image.getWidth());
    }

    // rescale the image the to next smaller power of 2
    if (VERBOSE_OUTPUT) {
      System.out.println("Rescale");
      System.out.flush();
    }
    int newSize = Integer.highestOneBit(size);
    float scaleFactor = (float)newSize / size;
    // System.out.printf("Rescale %d to %d: %f\n", size, newSize, scaleFactor);

    if (newSize != size) {
      AffineTransform rescaleXform = AffineTransform.getScaleInstance
        (scaleFactor, scaleFactor);
      AffineTransformOp rescaleOp = new AffineTransformOp
        (rescaleXform, AffineTransformOp.TYPE_BICUBIC);

      // must match the image type of the input image
      BufferedImage newImage = new BufferedImage
        (newSize, newSize, image.getType());

      rescaleOp.filter(image, newImage);
      image = newImage;
    }

    // convert to grayscale
    if (VERBOSE_OUTPUT) {
      System.out.println("Convert to grayscale");
      System.out.flush();
    }
    BufferedImage grayImage = new BufferedImage
      (newSize, newSize, BufferedImage.TYPE_BYTE_GRAY);
    Graphics g = grayImage.getGraphics();
    g.drawImage(image, 0, 0, null);
    g.dispose();
    image = grayImage;

    // enable this to output a grayscale jpg
    /*
    try {
      ImageIO.write(image, "jpeg", new File("output.jpg"));
    } catch (IOException e) {
      System.err.println("Error saving image: " + e);
    }
    */

    if (VERBOSE_OUTPUT) {
      System.out.println("Write data");
      System.out.flush();
    }
    try {
      long startNano = System.nanoTime();
      writeMatrix(image, outputFile, !textOutput);
      double elapsed = (System.nanoTime() - startNano) / 1e9;
      System.out.printf("%d x %d grayscale data written to \"%s\" " +
                        "in %.3f sec\n",
                        image.getWidth(), image.getHeight(), outputFile,
                        elapsed);
    } catch (IOException e) {
      System.err.println("Error writing to \"" + outputFile + "\": " + e);
    }
      
  }


  /* Read a text data file, convert to a grayscale JPEG. */
  public static void matrixToImage(String inputFile, String outputFile) {
    float data[];
    int width, height;

    try {
      Matrix rawData = readMatrix(inputFile);
      data = rawData.data;
      width = rawData.width;
      height = rawData.height;
    } catch (IOException e) {
      System.err.println("Error reading \"" + inputFile + "\": " + e);
      return;
    }

    // save the data as a grayscale image
    BufferedImage image = new BufferedImage
      (width, height, BufferedImage.TYPE_BYTE_GRAY);

    int i = 0;
    for (int y=0; y < height; y++) {
      for (int x=0; x < width; x++) {
        int color = (int)(data[i++] * 255 + .5);
        color = color | (color  << 8) | (color  << 16);
        image.setRGB(x, y, color);
      }
    }

    try {
      ImageIO.write(image, "jpeg", new File(outputFile));
    } catch (IOException e) {
      System.err.println("Error saving image: " + e);
    }

    System.out.printf("%d x %d image written to \"%s\".\n",
                      width, height, outputFile);
  }


  public static void writeMatrix(BufferedImage image, String outputFile,
                                  boolean binaryFormat)
    throws IOException {

    if (binaryFormat)
      writeMatrixBinary(image, outputFile);
    else
      writeMatrixText(image, outputFile);
  }


  public static Matrix readMatrix(String inputFile) throws IOException {
    if (isBinaryFile(inputFile)) {
      return readMatrixBinary(inputFile);
    } else {
      return readMatrixText(inputFile);
    }
  }    


  public static void writeMatrixText(BufferedImage image, String outputFile)
    throws IOException {

    DecimalFormat formatter = new DecimalFormat("0.00000");
    
    BufferedWriter out = new BufferedWriter
      (new FileWriter(outputFile));

    out.write(image.getWidth() + " cols " + image.getHeight() + " rows\n");

    for (int y=0; y < image.getHeight(); y++) {
      for (int x=0; x < image.getWidth(); x++) {
        if (x > 0) out.write(' ');
          
        int colorInt = image.getRGB(x, y) & 0xff;
        float colorFloat = colorInt / 255.0f;
        out.write(formatter.format(colorFloat));
      }
      out.write('\n');
    }
    out.close();
  }
      


  // Write in little-endian order
  private static void writeMatrixBinary(BufferedImage image, String outputFile)
    throws IOException {

    DataOutputStream out = new DataOutputStream
      (new BufferedOutputStream(new FileOutputStream(outputFile)));
    out.writeBytes("binary \n");

    out.writeInt(reverseBytes(image.getWidth()));
    out.writeInt(reverseBytes(image.getHeight()));

    for (int y=0; y < image.getHeight(); y++) {
      for (int x=0; x < image.getWidth(); x++) {
        int colorInt = image.getRGB(x, y) & 0xff;
        float colorFloat = colorInt / 255.0f;
        // float colorFloat = colorInt;

        out.writeInt(reverseBytes(Float.floatToRawIntBits(colorFloat)));

      }
    }
    out.close();
  }    


  private static final int reverseBytes(final int x) {
    return
      ((x & 0xff) << 24) |
      ((x & 0xff00) << 8) |
      ((x & 0xff0000) >>> 8) |
      ((x & 0xff000000) >>> 24);
  }


  public static boolean isBinaryFile(String filename) throws IOException {
    FileInputStream in = new FileInputStream(filename);
    byte buf[] = new byte[8];

    in.read(buf);
    in.close();

    String header = new String(buf);
    // System.out.println("header = \"" + header + "\"");
    return header.equals("binary \n");
  }


  private static Matrix readMatrixText(String inputFile) throws IOException {
    Scanner in = new Scanner(new File(inputFile));
    Matrix data = new Matrix();

    if (!in.hasNextInt()) {
      System.out.println("Bad header row");
      return null;
    }
    data.width = in.nextInt();

    if (!in.next().equals("cols")) {
      System.out.println("Bad header row");
      return null;
    }

    if (!in.hasNextInt()) {
      System.out.println("Bad header row");
      return null;
    }
    data.height = in.nextInt();


    if (!in.next().equals("rows")) {
      System.out.println("Bad header row");
      return null;
    }

    // read the float data
    int size = data.width * data.height;
    data.data = new float[size];

    for (int i=0; i < size; i++) {
      if (!in.hasNextFloat()) {
        System.out.printf("Error after %d values read.\n", i);
        return null;
      }
      data.data[i] = in.nextFloat();
    }

    in.close();

    return data;
  }


  private static Matrix readMatrixBinary(String inputFile) throws IOException {
    Matrix data = new Matrix();

    DataInputStream in = new DataInputStream
      (new BufferedInputStream(new FileInputStream(inputFile)));

    // skip over 8 bytes; should be "binary \n"
    in.readLong();

    data.width = reverseBytes(in.readInt());
    data.height = reverseBytes(in.readInt());
    if (VERBOSE_OUTPUT)
      System.out.printf("binary size %d x %d\n", data.width, data.height);
    int size = data.width * data.height;
    data.data = new float[size];
    for (int i=0; i < size; i++) {
      int tmp = in.readInt();
      data.data[i] = Float.intBitsToFloat(reverseBytes(tmp));
    }

    return data;
  }
}
    
