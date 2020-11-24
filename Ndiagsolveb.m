function [ x ] = Ndiagsolveb( cc, rr )
%NDIAGSOLVEB Solve a 2p+1-diagonal system
%   Solve a 2p+1-diagonal system on a bounded domain
%     Ax = rr
% The coefficients in row i are cc(1:2*p+1,i)

[Ndiag,N] = size(cc);
p = (Ndiag-1)/2;
d = p + 1;

% Eliminate terms below the diagonal
for k = 1:N-1
    % Using row k ...
    % ... eliminate from r further rows
    r = min(p,N-k);
    for i = 1:r
        % Eliminate from row kk
        kk = k + i;
        q = cc(d-i,kk)/cc(d,k);
        cc(d-i:d-i+r,kk) = cc(d-i:d-i+r,kk) - q*cc(d:d+r,k);
        rr(kk) = rr(kk) - q*rr(k);
    end
end

% Eliminate terms above the diagonal
for k = N:-1:2
    % Using row k ...
    % ...eliminate from r further rows
    r = min(p,k-1);
    for i = 1:r
        % Eliminate from row kk
        kk = k - i;
        q = cc(d+i,kk)/cc(d,k);
        cc(d+i,kk) = cc(d+i,kk) - q*cc(d,k);
        rr(kk) = rr(kk) - q*rr(k);
    end
end

% Now the problem is diagonal so solution is
x = rr./cc(d,:);


end