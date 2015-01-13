fid=fopen('Macular Cube 512x128.img', 'r');
%xsz=512;ysz=1024;zsz=128;
xsz=512;ysz=960;zsz=128;
Z=fread(fid,512*1024*128);tmp=reshape(Z,[512,1024,128]);
wo=tmp(:,65:1024,:);%Cut out portion with artifacts
  %Create structure with compression parameters
  S = s_init(8,.96,xsz,ysz,zsz,' VS_7.9');
  %Forward 3D Wavelet Transform
  w = comp(wo,S);
  %Threshold and quantize (loose information)
  [w,Xmax,th] = thresh(w,S.NZ,S.NX*S.NY,S.CUT);th
  w=reshape(w,[S.NZ,S.NX,S.NY]);
  wi=w(1:2:S.NZ,1:2:S.NX,1:2:S.NY);wi=wi(:);
  [partition,codebook]=lloyds(wi(:),256);
  wq=quantiz(w(:),partition,codebook);
  %Lossless compression (for testing only; will use bzip2 in C code)
  
%   [wqc,Res]=JPEGlike(0,wq);
%   Compression ratio (assuming input data are unsigned integers)
%   b2=size(wqc,1);b1=xsz*ysz*zsz;
%   ['Compression ratio = ',num2str(b1/b2)] 
%   ['Compressed file (Bytes) = ',num2str(b2)]
  %Dequantization
  w = codebook(wq+1);w=reshape(w,[S.NZ,S.NX*S.NY]);
  %Inverse 2D Wavelet Transform
  wr= ucomp(w,S); maxerror = max(abs(wo(:)-wr(:)));mse=sum((wo(:)-wr(:)).^2)/(xsz*ysz*zsz);
  PSNR=10*log10(255^2/mse);MW=max(wo(:));mW=min(wo(:));%Peak Signal-to-Noise Ratio
  for k=1:S.NZ
      ['slice = ',num2str(k)]
      ['Hit return for next slice']  
  figure(1);imagesc(wo(:,:,k));title('original');colormap(jet);
  figure(2);imagesc(wr(:,:,k),[mW,MW]);colormap('jet');title('reconstructed');
  pause
  end
  figure(3);H=vol3d('CData',wr,'texture','3D');colormap(gray(256));title('data cube');
  alphamap([0 linspace(0.1, 0, 255)]);
  axis equal off
  %set(gcf, 'color', 'w');
  view(3);[maxerror,mse,PSNR]
  % 2014 - S.E.Zarantonello