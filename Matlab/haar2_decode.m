function A = haar2_decode(B)

len = size(B,1);
lexp   = log2(len);

%build Haar filter matrix
T=[1  1; 1 -1];
I = eye(len);
H = kron(I(1:len/2,1:len/2),T)/sqrt(2);

%build permutation matrix
PT = I([1:2:len],:);PB = I([2:2:len],:);P=[PT;PB];
PD=P;

%decode Image
A=B;len2=1;
for j = 1:lexp
    len2 = 2*len2;
    AE  = A(1:len2,1:len2);
    PE  = [PT(1:len2/2,1:len2);PB(1:len2/2,1:len2)];
    HE  = H(1:len2,1:len2);
    AE  = HE'*PE'*AE*PE*HE;
    A(1:len2,1:len2)=AE;
endEnter file contents here
