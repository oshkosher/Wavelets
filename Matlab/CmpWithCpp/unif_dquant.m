function x = unif_dquant(y,amax,th,bits,thOnly)

% Inputs
%     y:   thresholded and quantized array
%  amax:   maximum absolute value of x
%    th:   threshold
%  bits:   bits/pixel for quantization
% Output
%     x:   dquantized array

base = 2^(bits-1)-1;
a    = abs(y);
[n1,n2] = size(y);
for j1=1:n1
for j2=1:n2
    if y(j1,j2)==0
      x(j1,j2)=0;
    else
        if thOnly ~= 1
            sign=y(j1,j2)/a(j1,j2);
            ax=a(j1,j2)*(amax)/base;
            x(j1,j2)=sign*ax;
        else
            x(j1,j2)=y(j1,j2);
        end
    end
end
end


