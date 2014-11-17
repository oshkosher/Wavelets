 function X = idata(FILE,slice,xsz,ysz)
 %FILE='C:\Work\DATA\TTest\TT.gt.15.H@';
 %xsz  = 201; ysz  = 101; slice = 1;
 %z=560; x=1120; y=450;
 fid=fopen(FILE,'r','l');
 w=(slice-1)*xsz*ysz*4; 
 fseek(fid,w,0);
 [X, count1] = fread(fid,[xsz,ysz],'real*4');
 %figure(slice);imagesc(X,[-.5,.5]);colormap hot(32)
 status=fclose(fid);
 A=imresize(X,[1024,1024],'bicubic');
 figure(slice);imagesc(X);colormap hot