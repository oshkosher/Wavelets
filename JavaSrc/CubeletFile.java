import java.nio.charset.Charset;
import java.io.*;
import com.google.protobuf.InvalidProtocolBufferException;

public class CubeletFile {

  private static final String HEADER_LINE = "SCU cubelets 1.0\n";
  private static final int HEADER_SIZE = 64;

  private static Charset ASCII;

  static {
    ASCII = Charset.forName("US-ASCII");
  }


  /** Encapsulate one cubelet.
   */
  public static class Cubelet {
    public int width, height, depth;
    public int xOffset, yOffset, zOffset;

    public enum DataType {UINT8, FLOAT32};
    public DataType datatype;

    public long dataFileOffset;

    public byte[] byteData;
    public float[] floatData;

    public void setSize(int w, int h, int d) {
      width = w;
      height = h;
      depth = d;
    }
    
    public void setOffset(int x, int y, int z) {
      xOffset = x;
      yOffset = y;
      zOffset = z;
    }

    public String toString() {
      return "Cubelet " + width + "x" + height + "x" + depth +
        ", offset " + xOffset + "," + yOffset + "," + zOffset;
    }
  }


  // test function - list the cubelets in a file
  public static void main2(String args[]) throws Exception {
    if (args.length != 1) printHelp();

    String filename = args[0];
    if (filename.equals("-")) filename = null;
    
    listCubelets(filename);
  }


  // test function - write simple cubelet file
  public static void main(String[] args) throws Exception {
    Cubelet cube = new Cubelet();
    float[] data = new float[100*100*100];

    String filename = args[0];
    if (filename.equals("-")) filename = null;

    Writer out = new Writer();
    out.open(filename);

    cube.datatype = Cubelet.DataType.FLOAT32;
    cube.floatData = data;

    cube.setSize(4, 5, 6);
    for (int i=0; i < 4*5*6; i++) data[i] = i*.1f;
    out.addCubelet(cube);
    
    cube.setSize(7, 8, 9);
    for (int i=0; i < 7*8*9; i++) data[i] = i*.2f;
    out.addCubelet(cube);

    out.close();
  }


  public static void listCubelets(String filename) throws Exception {
    Reader reader = new Reader();

    try {
      reader.open(filename);
    } catch (CubeletFormatException e) {
      System.out.println(e.getMessage());
      return;
    } catch (FileNotFoundException e) {
      System.out.println("\"" + filename + "\" not found.");
      return;
    } catch (IOException e) {
      System.out.println("Couldn't open cubelet file \"" + filename + 
                         "\":" + e);
      return;
    }

    Cubelet cube = new Cubelet();

    while (reader.next(cube)) {
      System.out.println(cube.toString());

      reader.getData(cube);

      if (cube.datatype == Cubelet.DataType.FLOAT32) {
        for (int y=0; y < cube.height; y++) {
          for (int x=0; x < cube.width; x++) {
            System.out.printf("%6.2f", cube.floatData[y*cube.width + x]);
          }
          System.out.println();
        }
      } else if (cube.datatype == Cubelet.DataType.UINT8) {
        for (int y=0; y < cube.height; y++) {
          for (int x=0; x < cube.width; x++) {
            System.out.printf("%4d", 0xff & cube.byteData[y*cube.width + x]);
          }
          System.out.println();
        }
      }
        
    }

    reader.close();
  }


  private static void printHelp() {
    System.out.println("\n  CubeletFile <file>\n" +
                       "  List the contents of the file.\n");
    System.exit(1);
  }
  
    
  // little endian
  public static long longFromByteArray(byte[] array, int offset) {
    long value = 0;
    for (int i=7; i >= 0; i--) {
      value = (value << 8) | (array[offset+i] & 0xff);
    }
    return value;
  }
    
  // little endian
  public static int intFromByteArray(byte[] array, int offset) {
    int value = 0;
    for (int i=3; i >= 0; i--) {
      value = (value << 8) | (array[offset+i] & 0xff);
    }
    return value;
  }
    
  // little endian
  public static void floatToByteArray(float value, byte[] array, int offset) {
    int intValue = Float.floatToIntBits(value);
    for (int i=0; i < 4; i++) {
      array[offset+i] = (byte) intValue;
      intValue >>>= 8;
    }
  }

  static class CubeletFormatException extends Exception {
    public CubeletFormatException(String filename) {
      super("\"" + filename + "\" is not a cubelet file.");
    }
  }


  static class CubeletStreamingException extends Exception {
    public CubeletStreamingException() {}
  }


  public static class Reader {
    private InputStream in;
    private long cubeletIndexOffset;
    private int dataSizeBytes = 0;
    private Cubelet.DataType currentCubeletDatatype;
    private boolean dataHasBeenRead = true;
    private boolean eofReached = false;

    // re-use a buffer for reading encoded messages
    // must be at least 4 bytes long, but will be grown as necessary
    private byte[] readBuffer = new byte[32];

    com.google.protobuf.Parser<WaveletCompress.CubeletBuffer> cubeletParser;


    public Reader() {
      cubeletParser = WaveletCompress.CubeletBuffer.getDefaultInstance()
        .getParserForType();
    }


    /** Opens the given file for reading.
        Throws CubeletFormatException if the file is not a cubelet file.
        Throws FileNotFoundException if the file is not found.
     */
    public void open(String filename)
      throws IOException, CubeletFormatException {

      if (in != null) in.close();
      
      if (filename == null) {
        in = System.in;
      } else {
        in = new BufferedInputStream
          (new FileInputStream(filename));
      }

      byte[] headerBuf = new byte[HEADER_SIZE];
      int bytesRead = in.read(headerBuf);
      if (bytesRead != HEADER_SIZE) {
        close();
        throw new IOException("Failed to read header");
      }

      String headerLine = new String(headerBuf, 0, HEADER_LINE.length(), ASCII);

      if (!headerLine.equals(HEADER_LINE)) {
        close();
        throw new CubeletFormatException(filename);
      }

      cubeletIndexOffset = longFromByteArray(headerBuf, 24);

      dataSizeBytes = 0;
      dataHasBeenRead = true;
      eofReached = false;
    }
    
    
    /** Reads a cubelet from the file, stores its metadata in 'cube',
        and returns true.
        Returns false when EOF reached or if the file has not been opened.
        Throws IOException on data format error or other IO error. */
    public boolean next(Cubelet cube) throws IOException {
      if (in == null || eofReached) return false;
      
      // if the data for this chunk hasn't been read, skip over it before
      // trying to read the next chunk.
      if (!dataHasBeenRead) {
        in.skip(dataSizeBytes);
      }

      WaveletCompress.CubeletBuffer cubeBuf = readCubeletBuffer();
      if (cubeBuf == null) {
        eofReached = true;
        return false;
      }

      cube.setSize(cubeBuf.getWidth(), cubeBuf.getHeight(), cubeBuf.getDepth());
      cube.setOffset(cubeBuf.getXOffset(), cubeBuf.getYOffset(),
                     cubeBuf.getZOffset());

      switch (cubeBuf.getDataType()) {
      case UINT8: cube.datatype = Cubelet.DataType.UINT8; break;
      case FLOAT32: cube.datatype = Cubelet.DataType.FLOAT32; break;
      }

      cube.dataFileOffset = 0;
      cube.byteData = null;
      cube.floatData = null;

      // save the size of the data for this buffer
      dataSizeBytes = cubeBuf.getByteCount();
      currentCubeletDatatype = cube.datatype;
      dataHasBeenRead = false;

      return true;
    }


    /** Get the data for the current cubelet.
        If the data is UINT8, store it in cube.byteData.
        If the data is FLOAT32, store it in cube.floatData.
        Throws CubeletStreamingException if the data for this cubelet has
        already been read (the input data may be streaming, so we can't
        back up to read the data again).
    */
    public void getData(Cubelet cube)
      throws IOException, CubeletStreamingException {

      if (dataSizeBytes == 0) return;
      
      if (dataHasBeenRead)
        throw new CubeletStreamingException();
      
      byte[] bytes = new byte[dataSizeBytes];
      int bytesRead = in.read(bytes);
      if (bytesRead != dataSizeBytes)
        throw new IOException("Failed to read cubelet data");

      dataHasBeenRead = true;
      
      if (currentCubeletDatatype == Cubelet.DataType.UINT8) {
        cube.byteData = bytes;
      }

      else if (currentCubeletDatatype == Cubelet.DataType.FLOAT32) {
        int count = dataSizeBytes / 4;
        float[] floats = new float[count];

        for (int i=0; i < count; i++) {
          // convert the bytes to ints, then to floats
          int tmp = intFromByteArray(bytes, i*4);
          floats[i] = Float.intBitsToFloat(tmp);
        }
        
        cube.floatData = floats;
      }

    }
    

    public void close() {
      if (in != null && in != System.in) {
        try {
          in.close();
        } catch (IOException e) {}
      }
      in = null;
    }


    // make sure the read buffer is at least this large
    private void readBufferReserve(int size) {
      if (readBuffer.length < size) {
        int newSize = Math.max(size, readBuffer.length*2);
        readBuffer = new byte[newSize];
      }
    }


    private WaveletCompress.CubeletBuffer readCubeletBuffer()
      throws IOException {

      // start by reading 4-byte length field
      int bytesRead = in.read(readBuffer, 0, 4);

      // end of file, just return null
      if (bytesRead != 4) return null;

      // decode the length
      int encodedLength = intFromByteArray(readBuffer, 0);

      // length of 0xFFFFFFFF marks EOF 
      if (encodedLength == -1) return null;
      
      readBufferReserve(encodedLength);

      // read the encoded data
      bytesRead = in.read(readBuffer, 0, encodedLength);
      if (bytesRead != encodedLength) return null;

      WaveletCompress.CubeletBuffer cubeletBuffer = null;

      // decode the metadata
      try {
        cubeletBuffer = cubeletParser.parseFrom(readBuffer, 0, encodedLength);
      } catch (InvalidProtocolBufferException e) {
        throw new IOException("Failed to decode metadata");
      }

      return cubeletBuffer;
    }
  }


  /**
     Pass-through OutputStream that counts the number of bytes
     written to it.
  */
  public static class CountingOutputStream extends OutputStream {
    private OutputStream out;
    private long byteCount;

    public long getByteCount() {
      return byteCount;
    }

    public CountingOutputStream(OutputStream out) {
      super();
      this.out = out;
    }

    public void write(byte[] b) throws IOException {
      out.write(b);
      byteCount += b.length;
    }

    public void write(byte[] b, int off, int len) throws IOException {
      out.write(b, off,len);
      byteCount += len;
    }
     
    public void write(int b) throws IOException {
      out.write(b);
      byteCount++;
    }
  }


  public static class Writer {

    private boolean isStdout;
    private CountingOutputStream out;
    
    public void open(String filename) throws IOException {
      if (out != null) close();
      
      if (filename != null) {
        isStdout = false;
        out = new CountingOutputStream
          (new BufferedOutputStream
           (new FileOutputStream(filename)));
      } else {
        out = new CountingOutputStream(System.out);
        isStdout = true;
      }
      
      byte[] headerBuf = new byte[HEADER_SIZE];

      // set the header line
      byte[] headerLine = HEADER_LINE.getBytes(ASCII);
      System.arraycopy(headerLine, 0, headerBuf, 0, headerLine.length);
      
      out.write(headerBuf);
    }

    public void addCubelet(Cubelet cube) throws IOException {
      if (out == null) throw new NullPointerException("File not open");

      if (cube.datatype == null)
        throw new NullPointerException("cubelet data type not set");

      int pixelSize = 0;
      int count = cube.width * cube.height * cube.depth;
      byte[] bytes = null;

      WaveletCompress.CubeletBuffer.DataType datatype = null;

      switch (cube.datatype) {
      case UINT8:
        datatype = WaveletCompress.CubeletBuffer.DataType.UINT8;
        if (cube.byteData == null)
          throw new NullPointerException("cubelet data not set");
        bytes = cube.byteData;
        pixelSize = 1;
        break;

      case FLOAT32:
        datatype = WaveletCompress.CubeletBuffer.DataType.FLOAT32;
        if (cube.floatData == null)
          throw new NullPointerException("cubelet data not set");

        pixelSize = 4;
        bytes = new byte[count * pixelSize];
        for (int i=0; i < count; i++)
          floatToByteArray(cube.floatData[i], bytes, i*4);
        break;

      }

      long cubeletFileOffset = out.getByteCount();

      WaveletCompress.CubeletBuffer cubeBuf
        = WaveletCompress.CubeletBuffer.newBuilder()
        .setWidth(cube.width)
        .setHeight(cube.height)
        .setDepth(cube.depth)
        .setXOffset(cube.xOffset)
        .setYOffset(cube.yOffset)
        .setYOffset(cube.zOffset)
        .setByteCount(pixelSize * count)
        .setDataFileOffset(cubeletFileOffset)
        .setDataType(datatype)
        .build();

      writeProtobuf(cubeBuf);
      
      // Set this field after writing out the metadata, because it's
      // only useful in the index at the end of the file.
      // cubeBuf.setDataFileOffset(cubeletFileOffset);
      
      // write the cubelet data
      out.write(bytes);
    }

    private void writeProtobuf(com.google.protobuf.GeneratedMessage m) {}

    public void close() throws IOException {
      if (out != null && !isStdout) {

        // XXX backpatch

        out.close();
      }
      out = null;
    }
  }
    
}

