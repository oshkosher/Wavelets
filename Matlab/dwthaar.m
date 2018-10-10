function r = dwthaar(Signal)

N = length(Signal);

s = zeros(1, N/2);
d = s;

for n=1:N/2
  s(n) = 1/2*(Signal(2*n-1) + Signal(2*n));
  d(n) = Signal(2*n-1) - s(n);
end

r = [s, d]
