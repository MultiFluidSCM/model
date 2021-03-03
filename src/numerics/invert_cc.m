function [ ] = invert_cc( cc, rr )
% Convert the compact representation cc into a full matrix
% and invert. Estimate which residuals contribute most to certain
% increments.
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

Ainv = inv(A);

disp('Row 16')
Ainv(16,1:27)
Ainv(16,1:27).*rr(1:27)
disp('Row 17')
Ainv(17,1:27)
Ainv(17,1:27).*rr(1:27)



end

