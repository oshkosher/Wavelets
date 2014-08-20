function th = threshold(A,cutoff)

[len,len] = size(A);

X  = sort(abs(A(:)));
th = X(floor(cutoff*len^2))+eps;
