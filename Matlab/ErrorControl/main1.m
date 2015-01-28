%Read input data
fid=fopen('Macular Cube 512x128.img', 'r');  %%%%IN PROGRESS
Z=fread(fid,512*1024*128,'uint8');tmp=reshape(Z,[512,1024,128]);
%Select portion of data
wo=tmp(:,65:1024,:); %Cut out portion of cube with artifacts
%Dimensions
xsz=512;ysz=960;zsz=128;m=xsz*ysz*zsz;

%Input parameters
ERRTYPE=' L1';max_rerr=0.15;WAVELET=' VS_7.9';
%Create structure
S = s_init(xsz,ysz,zsz,WAVELET);

%3D wavelet transform
w = fwt_3d(wo,S);

%%%%Error control loop begins
  %Goal: Highest compression ratio constrained by maximum relative error.
dcut=1;dbins=1;%Increments for variablse "cut" and "bins"
  %Norm of input data
nrm=get_nrm(wo,ERRTYPE);
% Step 1: Get initial values of "cut" and "bins"
   %Find maximum allowable cutoff without quantization
[max_cut,x]=get_max_cut(w,wo,nrm,max_rerr,ERRTYPE,S);
   %Find minimum number of bins consistent with max_cut-5*dcut
cut=max_cut-2*dcut;th=x(floor(cut/100*m));wt=thresh(w(:),m,th);
wt=reshape(wt,[zsz,xsz,ysz]);
min_bins=get_min_bins(wt,wo,nrm,max_rerr,ERRTYPE,S);
bins=min_bins+5*dbins;

% Step 2:"Hill-climbing" to get max "rat00" as function of "bins" and "cut" 
max_iter=10;           % maximum iterations
iter=1;rerr00=0;       % initial values for "while" loop
m=xsz*ysz*zsz;
while  iter < max_iter && rerr00 < max_rerr
  [rerr10,rat10]=rerr_ratio(cut+dcut,bins,wo,w,x,nrm,ERRTYPE,S);
  [rerr00,rat00]=rerr_ratio(cut,bins,wo,w,x,nrm,ERRTYPE,S);
  [rerr01,rat01]=rerr_ratio(cut,bins+dbins,wo,w,x,nrm,ERRTYPE,S);
  %Exit if error exceeds max rerr when increasing "cut" and decreasing "bins"
    if rerr01 > max_rerr && rerr10 > max_rerr
        break
    end
  %Climb to higher value of "rat00" 
    if rerr01 <= max_rerr && rerr10 <= max_rerr
        if rat10 >= rat00 && rat10 >= rat01
            cut=cut+dcut;rat00=rat10;
        elseif rat01 > rat00 && rat01 > rat10
            bins=bins-dbins;rat00=rat01;
        end
    elseif rerr01 > max_rerr && rerr10 <= max_rerr
        if rat10 >= rat00
            cut=cut+dcut;rat00=rat10;
        end
    elseif rerr01 <= max_rerr && rerr10 > max_rerr
        if rat01 >= rat00
            bins=bins-dbins;rat00=rat01;
        end
    [iter,cut,bins,rat00]
    iter=iter+1; 
    end
end 
%%%%Error control loop ends

%Apply optimal values of "cut" and "bins" for final result
th=x(floor(cut/100*m));[cut,bins,rat00]
wt=thresh(w(:),m,th);wi=reshape(wt,[zsz,xsz,ysz]);
wi=wi(1:4:zsz,1:4:xsz,1:4:ysz);
[partition,codebook]=lloyds(wi(:),bins);
wq=quantiz(wt(:),partition,codebook);wr=codebook(wq+1);
wr=reshape(wr,[zsz,xsz*ysz]);wr=iwt_3d(wr,S);
rerr00=get_rerr(wo,wr,nrm,ERRTYPE);   
[rat, original_bytes, compressed_bytes] = saveAndCompress(wq,'test');
rat00=m/compressed_bytes;
%Other outputs
maxerror = max(abs(uint8(wo(:)-wr(:))));mse=sum((wo(:)-wr(:)).^2)/m;
PSNR=10*log10(255^2/mse);ratio=m/compressed_bytes; 
'Peak Signal-to-Noise Ratio, comp ratio,rerrr'
[PSNR,ratio,rerr00]
for k=1:S.NZ
      ['slice = ',num2str(k)]
      ['Hit return for next slice']  
  figure(1);imagesc(wo(:,:,k));title('original');colormap(jet);
  figure(2);imagesc(uint8(wr(:,:,k)));colormap('jet');title('reconstructed');
  pause
end
% Using MSE=0.05: [PSNR,cut,bins,rerr00,rerr01,rerr10]
%                  28   92   97  0.0495 0.0495 0.0526
% 2015 - S.E.Zarantonello