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

    % java WaveletSampleImage Images/squirrel.jpg
     
     output:
       2048 x 2048 grayscale data written to "Images/squirrel.data" in 0.535 sec

   To convert the data file back to an image:

     % java WaveletSampleImage Images/squirrel.data squirrel2.jpg

     output:
       binary size 2048 x 2048
       2048 x 2048 image written to "squirrel2.jpg".

   The data files are in a binary format by default. Add the command
   line parameter "-text" to generate text files. See data_io.h for
   format details.
   
2) haar.cu

   This implements the Haar discrete wavelet transform on a 2-d data file.
   (or at least some broken version of the transform; it inverts correctly,
   but I don't know if it's doing the transform correctly).
   Give it an input data file and an output data file, and it performs
   the 2-d transform. Add the "-inverse" command line flag and it performs
   the inverse.
   
   It does every transform on the CPU and the GPU, and compares the
   timing and the results.  The implementations of each are in
   dwt_cpu.cc and dwt_gpu.cu.
   
   Example:

    % ./haar squirrel.data squirrel2.data
    Reading squirrel.data...2048 x 2048
    CPU: 130.714783 ms
    CUDA: 11.234272 ms
    Wrote squirrel2.data
   
    % ./haar -inverse squirrel2.data squirrel3.data
    Reading squirrel2.data...2048 x 2048
    CPU: 132.840942 ms
    CUDA: 11.967424 ms
    Wrote squirrel3.data
   
    % java WaveletSampleImage squirrel2.data
    binary size 2048 x 2048
    2048 x 2048 image written to "squirrel2.jpg".
    
    % java WaveletSampleImage squirrel3.data
    binary size 2048 x 2048
    2048 x 2048 image written to "squirrel3.jpg".

   
