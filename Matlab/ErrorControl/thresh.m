function y = thresh(w,m,th); %%%%IN PROGRESS
% w: input 1-D array 
% m: size of w
% th: threshold 
 
% Initialize output
x=abs(w);y=w; 
% Apply threshold 
for k = 1:m
  if(x(k)<= th) 
   y(k) = 0;
  end
end
%
% S.E.Zarantonello 2015