function [index,codebook,sign] = lyds_quant(B,th,bits)

% Inputs
%     x:   2-D array
%    th:   threshold
%  bits: number of bits for quantization
% Outputs:
%     y:    thresholded and quantized version of x
%     partition
%     codebook
[n1,n2]=size(B);base=2^(bits-1);q=zeros(n1,n2);sign=zeros(n1,n2);x=zeros(n1,n2);
for j1=1:n1
for j2=1:n2
    if abs(B(j1,j2))>th
        sign(j1,j2)=B(j1,j2)/abs(B(j1,j2));x(j1,j2)=abs(B(j1,j2));
    end
end
end
[partition,codebook]=lloyds(x(:),base);
index = quantiz(x(:),partition);index=reshape(index,n1,n2);