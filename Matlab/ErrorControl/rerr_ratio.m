function [rerr00,rat00]=rerr_ratio(cut,bins,wo,w,x,nrm,ERRTYPE,S)
m=length(wo(:));th=x(floor(cut/100*m));wt=thresh(w(:),m,th);
[xsz,ysz,zsz]=size(wo);
wi=reshape(wt,[zsz,xsz,ysz]);wi=wi(1:4:zsz,1:4:xsz,1:4:ysz);
[partition,codebook]=lloyds(wi(:),bins);
wq=quantiz(wt(:),partition,codebook);wr = codebook(wq+1);
wr=reshape(wr,[zsz,xsz*ysz]);wr=iwt_3d(wr,S);
rerr00=get_rerr(wo,wr,nrm,ERRTYPE);
[~,~, compressed_bytes] = saveAndCompress(wq,'test');
rat00=m/compressed_bytes;     
%2015 - S.E.Zarantonello