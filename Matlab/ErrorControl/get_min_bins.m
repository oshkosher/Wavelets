function min_bins=get_min_bins(wt,wo,nrm,max_rerr,ERRTYPE,S)
[zsz,xsz,ysz]=size(wt);wi=wt(1:4:zsz,1:4:xsz,1:4:ysz);
a=0;b=1024;min_bins=floor((b+a)/2);
for k=1:7
    [partition,codebook]=lloyds(wi(:),min_bins);
    wq=quantiz(wt(:),partition,codebook);wr = codebook(wq+1);
    wr=reshape(wr,[zsz,xsz*ysz]);wr=iwt_3d(wr,S);
    rerr=get_rerr(wo,wr,nrm,ERRTYPE);[rerr,min_bins]
    if rerr>=max_rerr
        a=min_bins;
    else
        b=min_bins;
    end
    if a==b
        min_bins=a;
        break
    end
    min_bins=floor((b+a)/2);
end
% S.E.Zarantonello 2015