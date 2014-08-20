function B=haar_encode(A);

len = size(A,1); % input image must be square matrix
lexp   = log2(len);
B = zeros(len,len);
%build Haar filter matrix
T=[1  1; 1 -1];
I = eye(len);
H = kron(I(1:len/2,1:len/2),T)/sqrt(2);

%build permutation matrix
PT = I([1:2:len],:);PB = I([2:2:len],:);P=[PT;PB];
PD=P;

%encode image
len2=len;
PD=P;B=A;
for j = 1:lexp
    PD  = [PT(1:len2/2,1:len2);PB(1:len2/2,1:len2)];
    BD  = B(1:len2,1:len2);
    HD  = H(1:len2,1:len2);
    BD  = PD*HD*BD*HD'*PD';
    B(1:len2,1:len2) = BD;
    len2 = len2/2;
end


