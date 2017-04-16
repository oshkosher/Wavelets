pkg load image                     % needed by Octave for rgb2gray

A=imread('Images/chicago_orig.jpg');  % convert jpeg image to RGB matrices
A=rgb2gray(A);                       % convert to grayscale matrix
dim1=256; dim2=256;                  % desired dimensions
A=imresize(A,[dim1,dim2],'bicubic'); % set A to desired dimensions
A=double(A);                         % convert A to floating points
figure(1);image(A);colormap(gray(256))% Input to Haar encode process
B=haar_encode(A);                    % One level Haar encoding
C=cudahaar(A);
figure(2);image(B);                  % Display result
colormap(gray(128));
figure(3);image(B(1:128,1:128)/2);   %
colormap(gray(256));
figure(4);image(C);                  % Display result
colormap(gray(128));
figure(5);image(C(1:128,1:128)/2);   %
colormap(gray(256));
