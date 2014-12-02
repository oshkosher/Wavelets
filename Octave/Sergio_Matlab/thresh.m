function [B,Xmax,Xmin] = thresh(A,len1,len2,cut);
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
    B=A;
    for i = 1:len1
    for j = 1:len2
        A1 = A(i,j); A2 = abs(A1);
        if(A2 <= Xmin) 
           B(i,j) = 0;
        end
    end
    end   
end
%
% Copyright (c) 2004. 3DGeo Development Inc.
% Comments ?  e-mail sergio@3dgeo.com