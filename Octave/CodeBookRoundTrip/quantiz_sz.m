function indx = quantiz_sz(w,p)
% Converts "w" into an array "indx" expressable with fewer bits 
% Inputs
%     w:   2-D array of wavelet coefficients
%     p:   partition p1 < p2 < . . . < pN    defining N+1 "bins"
% Outputs:
%     indx: quantized output (the array of bin numbers for "w")
% Note: 
%      1. "codebook" recreates lossy data: w = codebook(indx+1) 
%      2. there are N+1 bins, and N+1 corresponding codebook entries. 
%         bin(1)={w| w <= p(1)}, codebook(1)=value assigned to bin(1) 
%         . . .
%         bin(N)={w| p(N-1) < w <= p(N)}
%         bin(N+1)={w| p(N) < w}

[nRows, nCols] = size(w);   % to ensure proper output orientation
indx = zeros(nRows, nCols);
for i = 1 : length(p)
    indx = indx + (w > p(i));
end
