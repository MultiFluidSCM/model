function [ dBdz ] = upwind_bernoulli( dzp, w )
%   Calculate an upwind finite difference approximation to
%   d/dz (0.5*w^2)

nz = numel(dzp);
nzp = nz + 1;

% w squared
w2 = w.*w;

% Upwind index offset
ix = w < 0;
ix(1) = 1;
ix(end) = 0;


% Upwind gradient
dBdz = (w2([1:nzp] + ix) - w2([0:nz] + ix))./dzp([0:nz] + ix);


end

