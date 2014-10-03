function x = lyds_dquant(sign,codebook,index)

% Output:
%     x:    dquantized array
[n1,n2] = size(index); 
for j1=1:n1
for j2=1:n2
    x(j1,j2)=sign(j1,j2).*codebook(index(j1,j2)+1);
end
end