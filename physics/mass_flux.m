function [F, dFdma, dFdmb, dFdw] = mass_flux( grid, m, w )

% Compute mass flux given m and w
% Also terms needed for linearization

nz = grid.nz;
nzp = nz+1;

% For testing, use first order upwind
F(1) = 0;
dFdma(1) = 0;
dFdmb(1) = 0;
dFdw(1) = 0;
for k = 2:nz
    if w(k) > 0
        F(k) = w(k)*m(k-1);
        dFdma(k) = 0;
        dFdmb(k) = w(k);
        dFdw(k) = m(k-1);
    else
        F(k) = w(k)*m(k);
        dFdma(k) = w(k);
        dFdmb(k) = 0;
        dFdw(k) = m(k);
    end
end
F(nzp) = 0;
dFdma(nzp) = 0;
dFdmb(nzp) = 0;
dFdw(nzp) = 0;


end

