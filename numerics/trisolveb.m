function [ x ] = trisolveb( aa, bb, cc, r )
%TRISOLVEB Solve a tridiagonal system
%   Solve a tridiagonal system on a bounded domain
%     Ax = r
%   aa is the coefficient above the diagonal,
%   bb is the coefficient below the diagonal,
%   cc is the coefficient on the diagonal


nz = length(r);

% Forward elimination sweep
q(1) = -aa(1)/cc(1);
x(1) = r(1)/cc(1);
for k = 2:nz
  p = 1.0/(cc(k)+bb(k)*q(k-1));
  q(k) = -aa(k)*p;
  x(k) = (r(k)-bb(k)*x(k-1))*p;
end

% Backward pass
for k = nz-1:-1:1
  x(k) = x(k)+q(k)*x(k+1);
end

end