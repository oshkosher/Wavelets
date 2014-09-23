A=imread('museum.jpg');              % convert jpeg image to RGB matrices
A=rgb2gray(A);                       % convert to grayscale matrix
dim1=256; dim2=256;                  % desired dimensions
A=imresize(A,[dim1,dim2],'bicubic'); % set A to desired dimensions
A=double(A);                         % convert A to floating points
figure(1);image(A);colormap gray(256);%axis equal;axis([1 256 1 256]);

% Apply transformation
B = haar2_encode(A);
%B = d4_encode(A);

%figure(2);imagesc(B)

% Find threshold consistent with specified cutoff
th = threshold(B,cutoff);
% Create partition and codebook for uniform quantization
[partition,codebook] = unif_sz(B(:),th,bits);
% Quantization
indx = quantiz_sz(B,partition);

% Dequantization
BU = codebook(indx+1);b=find(BU==0);mb=length(b);
% Create 2D array
BU=reshape(BU,dim1,dim2);
% Inverse transformation
AU = haar2_decode(BU);
% AU = d4_decode(BU);
figure(3);image(AU);colormap gray(256);title('Uniform Quant');%axis equal;axis([1 256 1 256]);
        
% Calculate error and percent zeros
abs_err = sqrt(sum((A(:)-AU(:)).^2));normA = sqrt(sum(A(:).^2));
rel_err = abs_err/normA;
maxv=255;len=dim1;db=20*log10(maxv*len/abs_err);
pcnt_zeros =100*mb/len^2; pcnt_err = rel_err*100; [pcnt_zeros,pcnt_err,db]
