function [ z_tketop ] = TKEtop( zp, tke, thresh )
%TKEtop Determine the height at which the TKE first drops below
% some threshold

% Size of arrays
nz = numel(zp);

% Initialize search
% Make sure the final answer will not be below zp(1)
tkea = max(tke(1),thresh);

% Search
k = 1;
done = 0;
while k < nz & ~done
    k = k + 1;
    tkeb = tkea;
    tkea = tke(k);
    done = (tkea < thresh);
end

% Now interplate between last two values to
% get final answer
za = zp(k);
zb = zp(k-1);
a = (tkeb - thresh)/(tkeb - tkea);
b = 1 - a;
z_tketop = a*za + b*zb;

