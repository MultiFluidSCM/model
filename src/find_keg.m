function [ keg, dkegdwa, dkegdwc, dkegdwb ] = find_dkedz( dzp, w )
%   Calculate an upwind finite difference approximation to
%   the kinetic energy gradient
%   d/dz (0.5*w^2), along with terms needed in linearization

nz = numel(dzp);
nzp = nz + 1;

% KE = half w squared
w2 = 0.5*w.*w;

% Now compute its upwind gradient and derivatives
% wrt w at levels above, centred, and below
keg(1) = 0;
dkegdwa(1) = 0;
dkegdwc(1) = 0;
dkegdwb(1) = 0;
for k = 2:nz
    if w(k) > 0
        keg(k) = (w2(k  ) - w2(k-1))/dzp(k-1);
        dkegdwa(k) = 0;
        dkegdwc(k) =   w(k  )/dzp(k-1);
        dkegdwb(k) = - w(k-1)/dzp(k-1);
    else
        keg(k) = (w2(k+1) - w2(k  ))/dzp(k  );
        dkegdwa(k) =   w(k+1)/dzp(k  );
        dkegdwc(k) = - w(k  )/dzp(k  );
        dkegdwb(k) = 0;
    end
end
keg(nzp) = 0;
dkegdwa(nzp) = 0;
dkegdwc(nzp) = 0;
dkegdwb(nzp) = 0;
end

