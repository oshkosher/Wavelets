function x = log_dquant(y,lmax,th,bits)

% Inputs
%     y:    thresholded and quantized version of x
%  lmax:    maximum log2 of absolute value of x/th
%    th:    threshold
%  bits:    bits/pixel used for quantization
% Output:
%     x:    dquantized array

base    = 2^(bits-1)-1;
a       = abs(y);
[n1,n2] = size(y);
for j1=1:n1
for j2=1:n2
    if y(j1,j2)==0
      x(j1,j2)=0;
    else
      sign     =  y(j1,j2)/a(j1,j2);
      ln       =  a(j1,j2)*lmax/base;
      x(j1,j2) =  sign*th*2^ln;
    end
end
end
    

