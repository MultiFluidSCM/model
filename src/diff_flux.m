function [Dq, dDqdqa, dDqdqb, dDqdm ] = diff_flux( grid, kdiff, q , m, sfcflx)

% Compute tracer flux due to diffusion given tracer mixing
% ratio q, eddy diffusivity kdiff, and surface flux. Also
% compute terms needed for linearization

nz = grid.nz;
nzp = nz+1;

cc = m.*kdiff./grid.dzp;
% Diffusive flux
Dq = -cc.*(q(2:nzp) - q(1:nz));
% Derivatives wrt q above and below
dDqdqa = -cc;
dDqdqb =  cc;
% Derivative wrt m
dDqdm = - kdiff.*(q(2:nzp) - q(1:nz))./grid.dzp;

% Adjust first interior level to account for surface flux
% consistent with mass budget
% * This is now computed in find_surface_flux and added in tendencies
% Dq(1) = Dq(1) + grid.belowp(1)*sfcflx;

