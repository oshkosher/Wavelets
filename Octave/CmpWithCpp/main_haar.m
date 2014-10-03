cutoff = 0.93
bits = 8

%A=read_jpg('museum.jpg',256,256);len=256;
%A=read_jpg('amth_image.jpg',256,256);len=256;
%A=read_jpg('IMG1.jpg',256,256);len=256;
%A=read_bmp('test_image.bmp',256,256);len=256;
%A=imread('../../Images/museum.jpg','jpg');              % convert jpeg image to RGB matrices
%A=rgb2gray(A);                       % convert to grayscale matrix
%dim1=2048; dim2=2048;                  % desired dimensions
A = readDataFile('museum.data');
dim = size(A);
%A=imresize(A,[dim1,dim2],'bicubic'); % set A to desired dimensions
%A=double(A);                         % convert A to floating points

figNum = 2;
%figure(figNum);imagesc(A);figNum=figNum+1;colormap (gray(256));axis equal;axis([1 dim(1) 1 dim(2)]);

% Apply transformation
B = haar_encode(A);
%B = d4_encode(A);

%figure(figNum);imagesc(A);figNum=figNum+1;colormap (gray(256));axis equal;axis([1 dim(1) 1 dim(2)]);
Bmin = min(min(B))
Bmax = max(max(B))
figure(figNum);imagesc(B,[Bmin,Bmax]);title('encoded');figNum=figNum+1;colormap (gray(256));axis equal;axis([1 dim(1) 1 dim(2)]);

% Find threshold consistent with specified cutoff
th = threshold(B,cutoff)
[BL,lmax] = apply_thresh(B,th,bits);
%'Logarithmic quantization/dequantization'
%[BLQ,lmax] =  log_quant(B,th,bits);
%b=find(BLQ==0);mb=length(b);
%BL = log_dquant(BLQ,lmax,th,bits);

% Inverse transformation
AL = haar_decode(BL);
%AL = d4_decode(BL);
figure(figNum);image(AL);figNum=figNum+1;colormap (gray(256));title('Log Quant');axis equal;axis([1 dim(1) 1 dim(2)]);

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