fid=fopen('Macular Cube 512x128.img', 'r');
xsz=512;ysz=7*64;zsz=32;m=0;wo=zeros(xsz,ysz,zsz);
for k=1:32
    m=m+1;
    Z=fread(fid, [1024 1024], 'uint8');
    wo(:,:,m)=Z(1:512,65:512);
    %wo(:,:,65-m)=Z(1:512,513+64:1024);
end
  %Create structure with compression parameters
  S = s_init(8,.96,xsz,ysz,zsz,' VS_7.9');
  %Forward 3D Wavelet Transform
  w = comp(wo,S);
  %Threshold and quantize (loose information)
  [w,Xmax,th] = thresh(w,S.NZ,S.NX*S.NY,S.CUT);
  [wq,lmax] = log_quant(w,th,S.BITS);
  %Lossless compression (for testing only; will use bzip2 in C code)
  [wqc,Res]=JPEGlike(0,wq);
  %Compression ratio (assuming input data are floats)
  b2=size(wqc,1);b1=xsz*ysz*zsz;
  ['Compression ratio = ',num2str(b1/b2)] 
  ['Compressed file (Bytes) = ',num2str(b2)]
  %Dequantization
  %w = dquant(wq,S.NY_C,S.NX_C,S.BITS,Ymax,Ymin);
  w = log_dquant(wq,lmax,th,S.BITS);
  %Inverse 2D Wavelet Transform
  wr= ucomp(w,S); maxerror = max(abs(wo(:)-wr(:)))
  MW=max(wo(:));mW=min(wo(:));
  for k=1:S.NZ
      ['slice = ',num2str(k)]
      ['Hit return for next slice']  
  figure(1);imagesc(wo(:,:,k));title('original');colormap(jet);
  figure(2);imagesc(wr(:,:,k),[mW,MW]);colormap('jet');title('reconstructed');
  pause
  end
  figure(3);H=vol3d('CData',wo,'texture','3D');colormap(gray(256));title('data cube');
  alphamap([0 linspace(0.1, 0, 255)]);
  axis equal off
  %set(gcf, 'color', 'w');
  view(3);
  % 2014 - S.E.Zarantonello