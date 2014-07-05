Wavelets
========

SCU GPU Wavelet Project

To build:
  You'll need to have CUDA and Java installed.
  I've tested the tools on Linux and on Windows (inside Cygwin).
  At the command line, just run "make".
  
  I haven't tried integrating it in an IDE like Visual Studio yet.

There are two main tools:

1) WaveletSampleImage.java

   This converts between jpeg images and data files. The data files
   are two-dimensional arrays of floats, and they represent the image
   as a grayscale. Given an image file, it will crop it to be square,
   resize it down to the next smaller power of two, and convert it to
   grayscale.
   
   The first command line argument is the name of the input file.
   It expects a file with the suffix ".jpg", ".jpeg", or ".data".
   
   If you want to specify the name of the output file, make that the
   second command line argument. Otherwise, the suffix of the input file
   will be replaced automatically to name the output file.
   
   For example, to convert an image into a data file:

    % java WaveletSampleImage Images/haleakala_orig.jpg
    2048 x 2048 grayscale data written to "Images/haleakala_orig.data" in 0.444 sec

   To convert the data file back to an image:

    % java WaveletSampleImage Images/haleakala_orig.data haleakala2.jpg
    2048 x 2048 image written to "haleakala2.jpg".

   By default, the data files are in a binary format for best
   performance.  Add the command line parameter "-text" to generate
   text files. See data_io.h for format details.
   
2) haar.cu

   This implements the Haar discrete wavelet transform on a 2-d data file.
   Give it an input data file and an output data file, and it performs
   the 2-d transform. Add the "-inverse" command line flag and it performs
   the inverse. By default it will just do a single step of the transform.
   Add an integer argument after the filenames to specify a different number
   of steps.
   
   It does every transform on the CPU and the GPU, and compares the
   timing and the results.  The implementations of each are in
   dwt_cpu.cc and dwt_gpu.cu.
   
   Example:

    % java WaveletSampleImage Images/hubble4096.jpg hubble.data
    4096 x 4096 grayscale data written to "hubble.data" in 1.730 sec
    
    % ./haar hubble.data hubble2.data 3
    Reading hubble.data...4096 x 4096
    CPU: 539.509521 ms
    Times:
      Copy data to GPU:      10.391 ms
      Transform time:         7.597 ms (6 calls)
      Transpose time:         6.776 ms (6 calls)
      Copy data from GPU:    10.360 ms
    CUDA: 35.140831 ms
    Wrote hubble2.data
    
    % java WaveletSampleImage hubble2.data
    4096 x 4096 image written to "hubble2.jpg".
    
    % ./haar -inverse hubble2.data hubble3.data 3
    Reading hubble2.data...4096 x 4096
    CPU: 516.411377 ms
    Times:
      Copy data to GPU:      10.412 ms
      Transform time:         8.872 ms (6 calls)
      Transpose time:         7.022 ms (6 calls)
      Copy data from GPU:    10.511 ms
    CUDA: 36.833439 ms
    Wrote hubble3.data
    
    % java WaveletSampleImage hubble3.data
    4096 x 4096 image written to "hubble3.jpg".
   
