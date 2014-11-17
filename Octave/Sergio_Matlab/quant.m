function [B,Xmax,Xmin] = quant(A,len1,len2,nbits,cut);
% Thresholding & uniform quantization

% Initialize outputs
B = A; Xmin = 0; Xmax  = 0;
% Convert to linear array and sort abs. values
m = len1*len2; X = reshape(A,m,1); 
X = sort(abs(X)); Xmax = max(X);

if Xmax <= eps   % eps = machine "zero"       
    return    
else      
    Xmin = X(ceil(cut*m))+eps;
    clear X;
    % Threshold & quantize
    Base = 2^nbits-1;
    for i = 1:len1
    for j = 1:len2
        A1 = A(i,j); A2 = abs(A1);
        if(A2 <= Xmin) 
           B(i,j) = 0;
        else
           S = A1/A2;
           Q = round( Base*(A2-Xmin)/(Xmax-Xmin) );
           B(i,j) = S*Q;
        end
    end
    end   
end
% Copyright (c) 2014. S.E.Zarantonello