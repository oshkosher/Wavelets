function [y,lmax] = log_quant(x,th,bits)

% Inputs
%     x:   2-D array
%    th:   threshold
%  bits: number of bits for quantization
% Outputs:
%     y:    thresholded and quantized version of x
%  lmax:    maximum log2 of absolute value of x/th

base = 2^(bits-1)-1;
a    = abs(x);
%lmax = max(max(log2(a/th)));
lmax=log2(max(max(a))/th);
[n1,n2] = size(x);
for j1=1:n1
for j2=1:n2
    if a(j1,j2)<=th
      y(j1,j2)=0;
    else
      sign     =  x(j1,j2)/a(j1,j2);
      ln       =  log2(a(j1,j2)/th);
      q        =  ceil(base*ln/lmax);
      y(j1,j2) =  sign*q;
    end
end
end
