function [p,c] = unif_sz(w,th,bits)
%w=[-3,-1,-.1,0,.1,2,3,4,5,-2,-.2,2,0.03,0.04,0.5,9.01];th=1;bits=3;
%min_w<p(1)<p(2)...<p(N1)=-th<P(N1+1)=th<p(N1+2)<...<p(N1+N2)< max_w
%   1 bin   N1 bins                   N2 bins    1 bin
%   N = N1 + N2 +1 = M1 + M2 + 3
%indx = indx + (xorig > partition(i));%
% Inputs
%     w:   2-D array of wavelet coefficients
%     th:   threshold
% Outputs:
%     p=partition
%     c=codebook
% Note: this just build "c" and "p" for uniform quantization. 
% Call it first, then call quantiz_sz.m
N=2^bits;max_w=max(w(:));min_w=min(w(:));
if min_w < -th && max_w > th
    len1=-th-min_w; len2=max_w-th; len=len1+len2;
    N1=round((N-1)*len1/len);d1=len1/N1;N2=N-1-N1;d2=len2/N2;
    p(1)=min_w+d1;c(1)=min_w;
    for k=2:N1
        p(k)=min_w+k*d1;
        c(k)=(p(k)+p(k-1))/2;
    end
    p(N1+1)=th;codebook(N1+1)=0;
    for k=2:N2
        p(N1+k)=th+(k-1)*d2;
        c(N1+k)=(p(N1+k)+p(N1+k-1))/2;
    end
    c(N)=max_w;
end
