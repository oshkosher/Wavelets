import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import javax.swing.event.*;
import javax.swing.*;
import java.net.URL;
import java.io.IOException;
import java.io.File;
import java.io.FileNotFoundException;


/*******************************************************\ 
* CubeletViewer: application for viewing Cubelet files. *
\*******************************************************/
public class CubeletViewer extends JFrame {

  // make the font this many points larger than the default
  private static final float LABEL_FONT_SIZE_OFFSET = 3;

  // size to which the brightness icon will be scaled
  private static final int BRIGHTNESS_ICON_SIZE = 32;

  // store the brightness icon here 
  private static ImageIcon BRIGHTNESS_ICON;


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

    customizeUI();
    loadIcons();

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


  static void customizeUI() {
    Font labelFont = UIManager.getFont("Label.font");
    labelFont = labelFont.deriveFont
      ((float)(labelFont.getSize() + LABEL_FONT_SIZE_OFFSET));
    UIManager.put("Label.font", labelFont);
  }


  static void loadIcons() {
    URL iconURL = CubeletViewer.class.getResource("brightness.png");
    if (iconURL != null) {
      ImageIcon originalIcon = new ImageIcon(iconURL);
      Image tmpImage = originalIcon.getImage().getScaledInstance
        (BRIGHTNESS_ICON_SIZE, BRIGHTNESS_ICON_SIZE, Image.SCALE_SMOOTH);
      BRIGHTNESS_ICON = new ImageIcon(tmpImage);
    }
    // System.out.println("icon: " + brightnessIcon);
  }


  Cubelet cube;
  CubePanel cubePanel;
  JSlider zLevelSlider;
  JLabel zLevelLabel;
  JSlider brightnessSlider;

  public CubeletViewer(Cubelet cube, String filename) {
    super("Cubelet viewer");

    Dimension size;

    this.cube = cube;

    setFilenameTitle(filename);
    setKeyBindings();

    cubePanel = new CubePanel();
    cubePanel.setCube(cube);

    brightnessSlider = new JSlider(SwingConstants.VERTICAL, -255, 255, 0);
    // brightnessSlider.setBorder(BorderFactory.createEmptyBorder(5,5,0,10));
    size = brightnessSlider.getPreferredSize();
    size.height = 200;
    brightnessSlider.setPreferredSize(size);
    brightnessSlider.addChangeListener(new ChangeListener() {
        public void stateChanged(ChangeEvent ce) {
          setBrightness(brightnessSlider.getValue());
        }});
    JLabel brightnessLabel;
    if (BRIGHTNESS_ICON != null) {
      brightnessLabel = new JLabel(BRIGHTNESS_ICON);
    } else {
      brightnessLabel = new JLabel("Brt");
    }
    brightnessLabel.setAlignmentX(.5f);

    Box rightBox = Box.createVerticalBox();
    rightBox.add(brightnessSlider);
    rightBox.add(brightnessLabel);
    rightBox.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));


    zLevelSlider = new JSlider(0, cube.depth-1, 0);
    zLevelSlider.setPaintLabels(true);
    zLevelSlider.setPaintTicks(true);
    zLevelSlider.setMajorTickSpacing(50);
    zLevelSlider.addChangeListener(new ChangeListener() {
        public void stateChanged(ChangeEvent ce) {
          setZ(zLevelSlider.getValue());
        }});
    zLevelSlider.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    size = zLevelSlider.getPreferredSize();
    size.width = cube.width;
    zLevelSlider.setPreferredSize(size);

    zLevelLabel = new JLabel("9999", SwingConstants.RIGHT);
    zLevelLabel.setPreferredSize(zLevelLabel.getPreferredSize());
    zLevelLabel.setText("0");
    zLevelLabel.setAlignmentY(0);

    Box bottomBox = Box.createHorizontalBox();
    bottomBox.add(zLevelSlider);
    bottomBox.add(zLevelLabel);
    bottomBox.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));

    add(rightBox, "East");
    add(cubePanel, "Center");
    add(bottomBox, "South");

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
    zLevelLabel.setText(""+z);
  }


  public void setBrightness(int b) {
    cubePanel.setBrightness(b);
  }


  public void doExit() {
    System.exit(0);
  }

  public void setKeyBindings() {
    InputMap inputMap = getRootPane().
      getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
    ActionMap actionMap = getRootPane().getActionMap();

    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_Q, 0), "exit");
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0), "exit");

    actionMap.put("exit", new AbstractAction() {
	public void actionPerformed(ActionEvent e) {doExit();}});
  }


  /********************************************\ 
  * Point3D: encapsulate an (x,y,z) coordinate *
  \********************************************/
  public static class Point3D {
    public int x, y, z;

    public static Point3D parse(String s) {
      Point3D p = new Point3D();
      if (p.set(s))
        return p;
      else
        return null;
    }

    // parse a string in the form "1,2,3"
    // return true iff all three values are parsed successfully
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

  
  /*************************************************************\ 
  * CubePanel: component for displaying slices of a 3-d Cubelet *
  \*************************************************************/
  public static class CubePanel extends JPanel {
    private Cubelet cube;

    // linear remapping from whatever the data range is to 0..255
    private float dataScale = 1, dataOffset = 0;

    // z-level of the cube this is currently displayed
    private int z;

    // value added to each color component
    private int brightnessOffset = 0;

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

      if (cube.maxPossibleValue > 0) {
        dataOffset = 0;
        dataScale = 255.0f / cube.maxPossibleValue;
        return;
      }

      float min, max;
      if (cube.datatype == Cubelet.DataType.FLOAT32) {
        min = max = cube.floatData[offset];
        for (int i=offset; i < offset+length; i++) {
          float f = cube.floatData[i];
          if (f < min) {
            min = f;
          } else if (f > max) {
            max = f;
          }
        }
      } else {
        // cube.datatype == Cubelet.DataType.INT32
        min = max = cube.intData[offset];
        for (int i=offset; i < offset+length; i++) {
          float f = cube.intData[i];
          if (f < min) {
            min = f;
          } else if (f > max) {
            max = f;
          }
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

    public void setBrightness(int b) {
      if (b < -255) {
	brightnessOffset = -255;
      } else if (b > 255) {
	brightnessOffset = 255;
      } else {
	brightnessOffset = b;
      }

      refreshImage();
    }


    // Parameters have been updated--recompute the backing image in
    // a background thread
    private void refreshImage() {
      // image = createLevelImage(cube, z);
      // repaint();

      SwingWorker worker = new SwingWorker<BufferedImage, Void>() {
        public BufferedImage doInBackground() {
          return createLevelImage(cube, z, brightnessOffset);
        }

        public void done() {
          try {
            image = get();
            repaint();
          } catch (Exception e) {}
        }};
      worker.execute();
    }


    private BufferedImage createLevelImage(Cubelet cube, int z,
					   int brightnessOffset) {
      int dataIdx = z * cube.height * cube.width;

      rescaleData(dataIdx, cube.height * cube.width);

      BufferedImage img = new BufferedImage
        (cube.width, cube.height, BufferedImage.TYPE_BYTE_GRAY);

      int[] rgbArray = new int[cube.width * cube.height];
      
      if (cube.datatype == Cubelet.DataType.UINT8) {
        for (int i=0; i < rgbArray.length; i++) {
          int pixel = (cube.byteData[dataIdx+i] & 0xff) + brightnessOffset;
	  if (pixel < 0)
	    pixel = 0;
	  else if (pixel > 255)
	    pixel = 255;
          rgbArray[i] = pixel | (pixel << 8) | (pixel << 16);
        }
      }

      else if (cube.datatype == Cubelet.DataType.FLOAT32) {
        for (int i=0; i < rgbArray.length; i++) {
          float tmp = (cube.floatData[dataIdx+i] + dataOffset) * dataScale
	    + brightnessOffset;
          if (tmp < 0)
            tmp = 0;
          else if (tmp > 255)
            tmp = 255;
          int pixel = (int)tmp & 0xff;
          rgbArray[i] = pixel | (pixel << 8) | (pixel << 16);
        }
      }

      else if (cube.datatype == Cubelet.DataType.INT32) {
        for (int i=0; i < rgbArray.length; i++) {
          float tmp = (cube.intData[dataIdx+i] + dataOffset) * dataScale
	    + brightnessOffset;
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

      g.drawImage(image, 0, 0, null);
    }
  }

}
