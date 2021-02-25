function [Dk, dDkdka, dDkdkb ] = diff_flux_tke( grid, kdiff, tke , mbar)

% Compute vertical flux of tke due to diffusion given
% eddy diffusivity kdiff. Also compute terms needed for linearization.
% Zero flux BC at bottom and top.

nz = grid.nz;
nzp = nz+1;

cc = mbar.*kdiff./grid.dzw;

% Bottom level
Dk(1) = 0;
dDkdka(1) = 0;
dDkdkb(1) = 0;

% Diffusive flux
Dk(2:nz) = -cc(2:nz).*(tke(2:nz) - tke(1:nz-1));
% Derivatives wrt u above and below
dDkdka(2:nz) = -cc(2:nz);
dDkdkb(2:nz) =  cc(2:nz);

% Top level
Dk(nzp) = 0;
dDkdka(nzp) = 0;
dDkdkb(nzp) = 0;

