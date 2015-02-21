import java.nio.charset.Charset;
import java.io.*;
import java.util.ArrayList;
import com.google.protobuf.InvalidProtocolBufferException;

/**
   Java code for handling cubelet files.

   This uses the Google Protocol Buffer API, so when compiling
   or running it, you'll need "protobuf.jar" in the classpath.
   
   On Linux, if the libprotobuf-java package is installed,
   it will probably be in /usr/share/java/protobuf.jar.
   On Windows, a good place to store it would be
   at ./protobuf-2.6.0/protobuf.jar
*/
public class CubeletFile {

  private static final String HEADER_LINE = "SCU cubelets 1.0\n";
  private static final int HEADER_SIZE = 64;

  private static Charset ASCII;

  static {
    ASCII = Charset.forName("US-ASCII");
  }


  // test function - list the cubelets in a file
  public static void main(String args[]) throws Exception {
    if (args.length != 1) printHelp();

    String filename = args[0];
    if (filename.equals("-")) filename = null;
    
    listCubelets(filename);
  }


  // test function - write simple cubelet file
  public static void main2(String[] args) throws Exception {
    Cubelet cube = new Cubelet();
    float[] data = new float[100*100*100];

    String filename = args[0];
    if (filename.equals("-")) filename = null;

    Writer out = new Writer();
    out.open(filename);

    cube.datatype = Cubelet.DataType.FLOAT32;
    cube.floatData = data;

    cube.setSize(4, 5, 6);
    cube.setOffset(100, 101, 102);
    for (int i=0; i < 4*5*6; i++) data[i] = i*.25f;
    out.addCubelet(cube);
    
    cube.setSize(7, 8, 9);
    cube.setOffset(200, 201, 202);
    for (int i=0; i < 7*8*9; i++) data[i] = i*.5f;
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
    } catch (IOException e3) {
      System.out.println("Could not open cubelet file \"" + filename + 
                         "\":" + e3);
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
      } else if (cube.datatype == Cubelet.DataType.INT32) {
        for (int y=0; y < cube.height; y++) {
          for (int x=0; x < cube.width; x++) {
            System.out.printf("%7d", 0xff & cube.intData[y*cube.width + x]);
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
  
    
  /** Read a little endian long from a byte array. */
  public static long longFromByteArray(byte[] array, int offset) {
    long value = 0;
    for (int i=7; i >= 0; i--) {
      value = (value << 8) | (array[offset+i] & 0xff);
    }
    return value;
  }
    
  /** Read a little endian int from a byte array. */
  public static int intFromByteArray(final byte[] array, int offset) {
    int value = 0;
    for (int i=3; i >= 0; i--) {
      value = (value << 8) | (array[offset+i] & 0xff);
    }
    return value;
  }
    
  /** Write a little endian int to a byte array. */
  public static void intToByteArray(int value, byte[] array, int offset) {
    for (int i=0; i < 4; i++) {
      array[offset+i] = (byte) value;
      value >>>= 8;
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
      throws IOException, FileNotFoundException, CubeletFormatException {

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
        // skip() doesn't always skip the requested number of bytes,
        // so it must be called in a loop
        long bytesToSkip = dataSizeBytes;
        while (bytesToSkip > 0) {
          long bytesSkipped = in.skip(bytesToSkip);
          bytesToSkip -= bytesSkipped;
        }
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
      case INT32: cube.datatype = Cubelet.DataType.INT32; break;
      }

      switch (cubeBuf.getCompressionAlgorithm()) {
      case NONE: cube.compressionAlg = Cubelet.CompressionAlg.NONE; break;
      case ZLIB: cube.compressionAlg = Cubelet.CompressionAlg.ZLIB; break;
      case WAVELET: cube.compressionAlg = Cubelet.CompressionAlg.WAVELET; break;
      }

      cube.maxPossibleValue = cubeBuf.getMaximumValue();

      cube.dataFileOffset = 0;
      cube.byteData = null;
      cube.floatData = null;
      cube.intData = null;

      // save the size of the data for this buffer
      dataSizeBytes = cubeBuf.getByteCount();
      currentCubeletDatatype = cube.datatype;
      dataHasBeenRead = false;

      return true;
    }


    /** Get the data for the current cubelet.
        If the data is UINT8, store it in cube.byteData.
        If the data is FLOAT32, store it in cube.floatData.
        If the data is INT32, store it in cube.intData.
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

      else if (currentCubeletDatatype == Cubelet.DataType.INT32) {
        int count = dataSizeBytes / 4;
        int[] ints = new int[count];

        for (int i=0; i < count; i++) {
          // convert the bytes to ints, then to floats
          ints[i] = intFromByteArray(bytes, i*4);
        }
        
        cube.intData = ints;
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

    /** Read a CubeletBuffer from the input stream. */
    private WaveletCompress.CubeletBuffer readCubeletBuffer()
      throws IOException {

      // start by reading 4-byte length field
      int bytesRead = in.read(readBuffer, 0, 4);

      // end of file, just return null
      if (bytesRead != 4) return null;

      // decode the length
      int encodedLength = intFromByteArray(readBuffer, 0);
      /*
      System.out.printf("encodedLength %02x %02x %02x %02x = %d\n",
                        readBuffer[0]&0xff, readBuffer[1]&0xff,
                        readBuffer[2]&0xff, readBuffer[3]&0xff,
                        encodedLength);
      */

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

    public void flush() throws IOException {
      out.flush();
    }

    public void close() throws IOException {
      out.flush();
      out.close();
    }
  }


  /**
     Write a cubelet file. 
  */
  public static class Writer {

    private boolean isStdout;
    private CountingOutputStream out;
    ArrayList<WaveletCompress.CubeletBuffer> index
      = new ArrayList<WaveletCompress.CubeletBuffer>();

    public int width=0, height=0, depth=0;

    /** Open the given file for writing.
        If filename is null, use System.out.
        Write the header.
    */
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


    /**
       Write a cubelet (and its data) to the output stream.
    */
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
        for (int i=0; i < count; i++) {
          int tmp = Float.floatToIntBits(cube.floatData[i]);
          intToByteArray(tmp, bytes, i*4);
        }
        break;

      case INT32:
        datatype = WaveletCompress.CubeletBuffer.DataType.INT32;
        if (cube.intData == null)
          throw new NullPointerException("cubelet data not set");

        pixelSize = 4;
        bytes = new byte[count * pixelSize];
        for (int i=0; i < count; i++) {
          intToByteArray(cube.intData[i], bytes, i*4);
        }
        break;

      }

      WaveletCompress.CubeletBuffer cubeBuf
        = WaveletCompress.CubeletBuffer.newBuilder()
        .setWidth(cube.width)
        .setHeight(cube.height)
        .setDepth(cube.depth)
        .setXOffset(cube.xOffset)
        .setYOffset(cube.yOffset)
        .setZOffset(cube.zOffset)
        .setByteCount(pixelSize * count)
        .setDataType(datatype)
        .setMaximumValue(cube.maxPossibleValue)
        .build();

      // write the metadata
      writeProtobuf(cubeBuf);

      // save our current position in the output stream
      long cubeletFileOffset = out.getByteCount();

      cubeBuf = cubeBuf.toBuilder()
        .setDataFileOffset(cubeletFileOffset)
        .build();
      index.add(cubeBuf);
      // System.out.println("Write cubelet at offset " + cubeletFileOffset);

      // write the cubelet data
      out.write(bytes, 0, cubeBuf.getByteCount());
    }


    /**
       Writes a message to the output file, prefixed with its 4-byte length.
       Returns the number of bytes in the message.
    */
    private int writeProtobuf(com.google.protobuf.Message m)
      throws IOException {
      return writeProtobuf(m, true);
    }

                              
    /**
       Writes a message to the output file. If 'writeLength' is true,
       prefix the message with its 4-byte length.
       Returns the number of bytes in the message.
    */
    private int writeProtobuf(com.google.protobuf.Message m,
                              boolean writeLength)
      throws IOException {
      /*
        consider Parser<MessageType>.parseDelimitedFrom(InputStream),
        MessageLite.writeDelimitedTo(),
        CodedOutpuStream::writeVarint32ToArray
        (https://developers.google.com/protocol-buffers/docs/reference/cpp/google.protobuf.io.coded_stream#CodedOutputStream.WriteVarint32)
        and ReadVarint32FromArray (saved in cubelet_file.cc)
      */

      // get the encoded length of the message
      int messageLen = m.getSerializedSize();
      
      // write the length as a 4-byte int
      if (writeLength) {
        byte[] lengthBuf = new byte[4];
        intToByteArray(messageLen, lengthBuf, 0);
        out.write(lengthBuf);
      }

      // write the encoded message
      m.writeTo(out);

      return messageLen;
    }


    /**
       Write the footer for the cubelet file and close the file.
    */
    public void close() throws IOException {
      if (out == null) return;

      // mark the end of cubelets with a length of 0xFFFFFFFF
      byte[] footerBuf = new byte[8];
      intToByteArray(-1, footerBuf, 0);
      out.write(footerBuf, 0, 4);

      // build the footer
      WaveletCompress.CubeletIndexBuffer.Builder indexBuilder
        = WaveletCompress.CubeletIndexBuffer.newBuilder();
      for (int i=0; i < index.size(); i++)
        indexBuilder.addCubelets(index.get(i));

      if (width > 0)  indexBuilder.setWidth(width);
      if (height > 0) indexBuilder.setHeight(height);
      if (depth > 0)  indexBuilder.setDepth(depth);
        
      WaveletCompress.CubeletIndexBuffer indexBuf = indexBuilder.build();
        
      // write the footer
      int indexLen = writeProtobuf(indexBuf, false);

      // at the end of the file, write an offset to the location of the
      // footer data, as well as a special byte sequence so we can
      // distinguish between a completed file and one that got truncated.
      intToByteArray(indexLen, footerBuf, 0);
      footerBuf[4] = 'c';
      footerBuf[5] = 'u';
      footerBuf[6] = 'b';
      footerBuf[7] = 'e';
      out.write(footerBuf);

      out.flush();

      if (!isStdout) out.close();

    }
  }
    
}

