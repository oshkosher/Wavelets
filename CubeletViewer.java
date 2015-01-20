import java.awt.*;
import java.awt.image.*;
import javax.swing.event.*;
import javax.swing.*;
import java.io.IOException;
import java.io.File;
import java.io.FileNotFoundException;

public class CubeletViewer extends JFrame {

  static class Point3D {
    public int x, y, z;

    public static Point3D parse(String s) {
      Point3D p = new Point3D();
      if (p.set(s))
        return p;
      else
        return null;
    }
    
    public boolean set(String s) {
      int start = 0, end = 0;

      try {
        end = s.indexOf(',');
        if (end == -1) return false;
        x = Integer.parseInt(s.substring(start, end));
        start = end+1;
        end = s.indexOf(',', start);
        if (end == -1) return false;
        y = Integer.parseInt(s.substring(start, end));
        z = Integer.parseInt(s.substring(end+1));
        return true;
      } catch (NumberFormatException e) {
        return false;
      }
    }

    public void set(int x, int y, int z) {
      this.x = x; this.y = y; this.z = z;
    }

    public boolean equals(int x, int y, int z) {
      return this.x == x && this.y == y && this.z == z;
    }

    public String toString() {
      return x + "," + y + "," + z;
    }
  }

  public static void main(String[] args) {

    if (args.length < 1 || args.length > 2) printHelp();

    String cubeFile = args[0];
    Point3D cubeId = null;

    if (args.length > 1) {
      cubeId = Point3D.parse(args[1]);
      if (cubeId == null) {
        System.out.println("Invalid cublet index: \"" + args[1] + "\"");
        return;
      }
    }
    
    CubeletFile.Reader reader = new CubeletFile.Reader();
    try {
      reader.open(cubeFile);
    } catch (FileNotFoundException e) {
      System.out.println("\"" + cubeFile + "\" not found");
                         
      return;
    } catch (Exception e) {
      System.out.println("Couldn't open cubelet file \"" + cubeFile + 
                         "\":" + e);
      return;
    }

    Cubelet cube = new Cubelet();

    while (true) {
      boolean hasNext = false;

      try {
        hasNext = reader.next(cube);
      } catch (IOException e) {
        System.out.println("Error reading cubelet file: " + e);
        return;
      }

      if (!hasNext) {
        cube = null;
        break;
      }

      // only uncompressed 8-bit int or 32-bit float supported
      if (cube.datatype != Cubelet.DataType.UINT8 &&
          cube.datatype != Cubelet.DataType.FLOAT32) continue;

      if (cube.compressionAlg != Cubelet.CompressionAlg.NONE) continue;

      // no ID was specified; use the first cubelet
      if (cubeId == null) break;

      if (cube.xOffset == cubeId.x &&
          cube.yOffset == cubeId.y &&
          cube.zOffset == cubeId.z) break;
    }

    if (cube == null) {
      System.out.println("No matching cubelet found.");
      return;
    }

    try {
      reader.getData(cube);
    } catch (Exception e) {
      System.out.println("Failed to read data for cubelet " + cubeId
                         + ": " + e);
      return;
    }

    /*
    if (cube.byteData.length != cube.width * cube.height * cube.depth) {
      System.out.println("Cubelet has " + cube.byteData.length + " bytes of " +
                         "data, expected " +
                         (cube.width * cube.height * cube.depth));
      return;
    }
    */

    CubeletViewer viewer = new CubeletViewer(cube, cubeFile);
    viewer.setVisible(true);

  }


  public static void printHelp() {
    System.out.println
      ("\n  CubeletViewer <cubelet file> [<x>,<y>,<z>]\n" +
       "  For example:\n" +
       "    CubeletViewer bigdatafile.cube 0,256,1024\n" +
       "  If the cubelet ID is omitted, the first cubelet in the file is used.\n");
    System.exit(1);
  }


  Cubelet cube;
  CubePanel cubePanel;
  JSlider zLevelSlider;

  public CubeletViewer(Cubelet cube, String filename) {
    super("Cubelet viewer");
    this.cube = cube;

    setFilenameTitle(filename);

    cubePanel = new CubePanel();
    cubePanel.setCube(cube);

    zLevelSlider = new JSlider(0, cube.depth-1, 0);
    zLevelSlider.setPaintLabels(true);
    zLevelSlider.setPaintTicks(true);
    zLevelSlider.setMajorTickSpacing(25);
    zLevelSlider.addChangeListener(new ChangeListener() {
        public void stateChanged(ChangeEvent ce) {
          setZ(zLevelSlider.getValue());
        }});
    Dimension size = zLevelSlider.getPreferredSize();
    size.width = cube.width;
    zLevelSlider.setPreferredSize(size);

    add(cubePanel, "Center");
    add(zLevelSlider, "South");

    pack();
    setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
    // setSize(new Dimension(500, 500));
    setLocationRelativeTo(null);
  }

  void setFilenameTitle(String filename) {
    File f = new File(filename);
    setTitle(f.getName());
  }

  public void setZ(int z) {
    cubePanel.setZ(z);
  }

  
  public static class CubePanel extends JPanel {
    private Cubelet cube;

    // linear remapping from whatever the data range is to 0..255
    private float dataScale = 1, dataOffset = 0;

    // z-level of the cube this is currently displayed
    private int z;

    // image created from the current z-level of the cube;
    BufferedImage image;

    public CubePanel() {
      super();
    }

    public void setCube(Cubelet cube) {
      this.cube = cube;
      z = 0;
      setPreferredSize(new Dimension(cube.width, cube.height));
      // if the data is not 8-bit integers, set constants to rescale it
      refreshImage();
    }

    void rescaleData(int offset, int length) {
      if (cube.datatype == Cubelet.DataType.UINT8) return;

      float min, max;
      min = max = cube.floatData[offset];
      for (int i=offset; i < offset+length; i++) {
        float f = cube.floatData[i];
        if (f < min) {
          min = f;
        } else if (f > max) {
          max = f;
        }
      }

      // System.out.println("min = " + min + ", max = " + max);

      if (min == max) {
        dataScale = 0;
        dataOffset = 0;
      } else {
        dataScale = 255.0f / (max - min);
        dataOffset = -min;
      }
    }

    
    public void setZ(int z) {
      if (z < 0) {
        this.z = 0;
      } else if (z >= cube.depth) {
        this.z = cube.depth-1;
      } else {
        this.z = z;
      }

      refreshImage();
    }

    private void refreshImage() {
      // image = createLevelImage(cube, z);
      // repaint();

      SwingWorker worker = new SwingWorker<BufferedImage, Void>() {
        public BufferedImage doInBackground() {
          return createLevelImage(cube, z);
        }

        public void done() {
          try {
            image = get();
            repaint();
          } catch (Exception e) {}
        }};
      worker.execute();
    }
    
    private BufferedImage createLevelImage(Cubelet cube, int z) {
      int dataIdx = z * cube.height * cube.width;

      rescaleData(dataIdx, cube.height * cube.width);

      BufferedImage img = new BufferedImage
        (cube.width, cube.height, BufferedImage.TYPE_BYTE_GRAY);

      int[] rgbArray = new int[cube.width * cube.height];
      
      if (cube.datatype == Cubelet.DataType.UINT8) {
        for (int i=0; i < rgbArray.length; i++) {
          int pixel = cube.byteData[dataIdx+i] & 0xff;
          rgbArray[i] = pixel | (pixel << 8) | (pixel << 16);
        }
      }

      else if (cube.datatype == Cubelet.DataType.FLOAT32) {
        for (int i=0; i < rgbArray.length; i++) {
          float tmp = (cube.floatData[dataIdx+i] + dataOffset) * dataScale;
          if (tmp < 0)
            tmp = 0;
          else if (tmp > 255)
            tmp = 255;
          int pixel = (int)tmp & 0xff;
          rgbArray[i] = pixel | (pixel << 8) | (pixel << 16);
        }
      }

      else {
        System.err.println("Unknown cubelet data type: " + cube.datatype);
        throw new NullPointerException();
      }
        

      img.setRGB(0, 0, cube.width, cube.height, rgbArray, 0, cube.width);
        
      /*
      for (int y=0; y < cube.height; y++) {
        for (int x=0; x < cube.width; x++) {

          // strip off high bits to cast -127 into 255
          int pixel = cube.byteData[dataIdx++] & 0xff;
          int rgb = pixel | (pixel << 8) | (pixel << 16);

          // System.out.println(pixel + " -> " + rgb);

          img.setRGB(x, y, rgb);

          // make a checkerboard, just for testing
        }
      }
      */

          /*
          int bit1 = y & (1<<5);
          int bit2 = x & (1<<5);
          if ((bit1 ^ bit2) == 0) {
            img.setRGB(x, y, 0x00000000);
          } else {
            img.setRGB(x, y, 0x00ffffff);
          }
          */

      return img;
    }

    public void paintComponent(Graphics g1) {
      super.paintComponent(g1);

      Graphics2D g = (Graphics2D) g1;

      // scale to fit the smallest dimension
      Dimension size = getSize();
      if (size.width <= 0 || size.height <= 0) return;
      double xScale = (double)size.width / cube.width;
      double yScale = (double)size.height / cube.height;
      double scale = Math.min(xScale, yScale);
      g.scale(scale, scale);

      // System.out.println("draw image " + image);
      g.drawImage(image, 0, 0, null);

      // g.setColor(Color.BLUE);
      // g.fillRect(10, 10, getSize().width-20, getSize().height-20);
    }
  }

}
