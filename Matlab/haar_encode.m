function B=haar_encode(A);
len = size(A,1); B = zeros(len,len);
%build Haar filter matrix
T=[1  1; 1 -1]; I = eye(len);
H = kron(I(1:len/2,1:len/2),T);
%build permutation matrix
PT = I([1:2:len],:);PB = I([2:2:len],:);P=[PT;PB];
%encode image one level
B=P*H*A*H'*P';B=B/2;
  



