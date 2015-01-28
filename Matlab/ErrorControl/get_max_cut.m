function [max_cut,x]=get_max_cut(w,wo,nrm,max_rerr,ERRTYPE,S) 
x=sort(abs(w(:)));m=length(x);
[xsz,ysz,zsz]=size(wo);m=length(wo(:));
a=0;b=100;
for k=1:7
    max_cut=floor((b+a)/2);
    th=x(ceil(max_cut/100*m));
    wt=thresh(w(:),m,th);wt=reshape(wt,[zsz,xsz*ysz]);ww=iwt_3d(wt,S);
    rerr=get_rerr(wo,ww,nrm,ERRTYPE);[rerr,max_cut]
    if rerr >=max_rerr
        b=max_cut;
    else
        a=max_cut;
    end
    if b-a <=2
        break
    end
end
max_cut=a
% S.E.Zarantonello 2015