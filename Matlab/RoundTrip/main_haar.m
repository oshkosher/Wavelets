%A=read_jpg('museum.jpg',256,256);len=256;
%A=read_jpg('amth_image.jpg',256,256);len=256;
%A=read_jpg('IMG1.jpg',256,256);len=256;
%A=read_bmp('test_image.bmp',256,256);len=256;
A=imread('../../Images/museum.jpg','jpg');              % convert jpeg image to RGB matrices
A=rgb2gray(A);                       % convert to grayscale matrix
dim1=256; dim2=256;                  % desired dimensions
A=imresize(A,[dim1,dim2],'bicubic'); % set A to desired dimensions
A=double(A);                         % convert A to floating points
figure(1);image(A);colormap gray(256);%axis equal;axis([1 256 1 256]);

% Apply transformation
B = haar_encode(A);
%B = d4_encode(A);

%figure(2);imagesc(B)

% Find threshold consistent with specified cutoff
th = threshold(B,cutoff);

'Logarithmic quantization/dequantization'
[BLQ,lmax] =  log_quant(B,th,bits);
b=find(BLQ==0);mb=length(b);
BL = log_dquant(BLQ,lmax,th,bits);

% Inverse transformation
AL = haar_decode(BL);
%AL = d4_decode(BL);
figure(2);image(AL);colormap gray(256);title('Log Quant');%axis equal;axis([1 256 1 256]);

% Calculate error and percent zeros
% abs_err = sum((A(:)-AL(:))^2;energy =sum(A(:).^2);

abs_err = sqrt(sum((A(:)-AL(:)).^2));
normA = sqrt(sum(A(:).^2));
rel_err = abs_err/normA;len=dim1;
maxv=255;db=20*log10(maxv*len/abs_err);PSNR=db
pcnt_zeros =100*mb/len^2; pcnt_err = rel_err*100; [pcnt_zeros,pcnt_err,PSNR]

%'Uniform quantization/dequantization'
% [BUQ,bmax] = unif_quant(B,th,bits);
% b=find(BUQ==0);mb=length(b);
% BU = unif_dquant(BUQ,bmax,th,bits);

% Inverse transformation
%AU = haar_decode(BU);
% AU = d4_decode(BU);
% figure(2);image(AU);colormap gray(256);title('Uniform Quant');%axis equal;axis([1 256 1 256]);
        
% Calculate error and percent zeros
% abs_err = sqrt(sum((A(:)-AU(:)).^2));
% rel_err = abs_err/normA;
% db=20*log10(maxv*len/abs_err);
% pcnt_zeros =100*mb/len^2; pcnt_err = rel_err*100; [pcnt_zeros,pcnt_err,db]