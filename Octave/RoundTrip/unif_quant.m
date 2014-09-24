function [y,amax] = unif_quant(x,th,bits)

% Inputs
%     x:   2-D array of wavelet coefficients
%    th:   threshold
%  bits: number of bits for quantization
% Outputs:
%     y:    thresholded and quantized version of x
%  amax:    maximum absolute value of x

base = 2^(bits-1)-1;
a    = abs(x);
amax = max(max(a));
[n1,n2] = size(x);
for j1=1:n1
for j2=1:n2
    if a(j1,j2)<=th
      y(j1,j2)=0;
    else
      sign=x(j1,j2)/a(j1,j2);
      y(j1,j2)=sign*ceil( base*(a(j1,j2)-th)/(amax-th) );
    end
end
end


