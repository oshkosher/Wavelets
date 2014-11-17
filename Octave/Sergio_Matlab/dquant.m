function B = dquant(A,len1,len2,nbits,Xmax,Xmin);

% Convert back from integers to floating point numbers

B = zeros(len1,len2); 
% Base = 2^(nbits-1);
Base = 2^nbits-1;
for i = 1:len1
for j = 1:len2
   A1 = A(i,j); A2 = abs(A1);
   if A2 == 0  
            B(i,j) = 0;
   else
           S = A1/A2; 
           Q = A2/Base*(Xmax-Xmin) + Xmin;%Xmax,Xmin,Base
           B(i,j)= S*Q;
   end
end
end
%Q = round( Base*(A2-Xmin)/(Xmax-Xmin) )
% Copyright (c) 2014. S.E.Zarantonello