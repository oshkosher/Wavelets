
PREREQUISITES:
---------------

- CUDA
- Octave Image package.


FILES:
-------

- cudahaar.oct:   An Octave file (like .m for matlab) that implements haar wavelets using CUDA.
- haar_encode.m:  Implements haar wavelets only using Octave.
- libhaar.so:     Shared library with the implementation in CUDA.
- main.m:         Example script to compare the haar_encode and cudahaar.
- README.txt:     this file


HOW TO USE:
------------

Fisrt, Octave must be able to find libhaar.so, this can be done in many ways but the easiest probably is to copy libhaar.so to /usr/lib.

Then executing inside Octave main.m should display 5 images:

	figure 1 is the original image in gray scale and resized.
	figure 2 and 3 are generated with the output of haar_encode.
	figure 4 and 5 are generated with the output of cudahaar.

The cudahaar functions is as follows:

	Output = cudahaar(Input[,step,inverted,text,blockSize]);

whith:

	Input:  A squared float point matrix.
	Output: A float point matrix with the same size as Input.
	step:     Optional parameter with the number of encoding levels. Default 1.
	inverted: A boolean value indicating false for encoding, true for decoding. Defaul false.
	blockSize: A integer indicating the size of the CUDA thread block. Default -1 that means that will be calculated as a function of the size of the input. 

