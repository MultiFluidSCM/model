function [ r ] = Ndiagmult( cc, x )
% Multiply x an N-diagonal matrix, represented compactly
% by the coefficients cc
% The coefficients in row i are cc(1:2*p+1,i)

[Ndiag,N] = size(cc);
p = (Ndiag-1)/2;
d = p + 1;

for k = 1:N
    dl = max(-p,1-k);
    dr = min( p,N-k);
    r(k) = x(k+dl:k+dr)*cc(d+dl:d+dr,k);
end

