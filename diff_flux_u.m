function [Du, dDudua, dDudub ] = diff_flux_u( grid, kdiff, u , mbar, speed, k_vonkarman, zrough)

% Compute vertical flux of u and v due to diffusion given
% eddy diffusivity kdiff. Also compute terms needed for linearization.
% A no slip BC is assumed at the bottom and a free slip BC at the top

nz = grid.nz;
nzp = nz+1;

cc = mbar.*kdiff./grid.dzw;

% Bottom level
coeff = mbar(1)*(k_vonkarman/log((grid.zp(1) + zrough)/zrough))^2;
Du(1) = -coeff*speed*u(1);
if speed == 0
    dDudua(1) = 0;
else
    dDudua(1) = -coeff*(speed + (u(1)^2)/speed);
end
dDudub(1) =  0;

% Diffusive flux
Du(2:nz) = -cc(2:nz).*(u(2:nz) - u(1:nz-1));
% Derivatives wrt u above and below
dDudua(2:nz) = -cc(2:nz);
dDudub(2:nz) =  cc(2:nz);

% Top level
Du(nzp) = 0;
dDudua(nzp) = 0;
dDudub(nzp) = 0;

