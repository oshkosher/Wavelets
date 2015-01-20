/** Encapsulate one cubelet.
 */
public class Cubelet {
  public int width, height, depth;
  public int xOffset, yOffset, zOffset;

  public enum DataType {UINT8, FLOAT32};
  public DataType datatype;

  public long dataFileOffset;

  public enum CompressionAlg {NONE, ZLIB, WAVELET};
  public CompressionAlg compressionAlg;

  // These correspond to the setting of 'datatype'.
  // For example, if datatype==UINT8, the data is assumed
  // to be in 'byteData', and 'floatData' is ignored.
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
