function [ x ] = Ndiagsolvex( cc, rr )
%Ndiagsolvex Solve an N-diagonal linear system
%   This version provides a sanity check on the version Ndiagsolveb
%   by converting the compact representation of the matrix to a full
%   explicit representation and using Matlabs inbuilt linear solver
% The coefficients in row i are cc(1:2*p+1,i)

[Ndiag,N] = size(cc);
p = (Ndiag-1)/2;
d = p + 1;

% Construct explicit matrix
A = zeros(N,N);
for i = 1:N
    j1 = max(d+1-i,1);
    j2 = min(d+N-i,Ndiag);
    k1 = max(i-p,1);
    k2 = min(i+p,N);
    A(i,k1:k2) = cc(j1:j2,i);
end
disp('Condition number of matrix:')
cond(A)

x = A\rr';
x = x';



end

