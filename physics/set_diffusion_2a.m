function [kdifft1,kdiffq1,kdiffw1,kdifft2,kdiffq2,kdiffw2,kdiffu1,kdiffu2,kdifftke1,kdifftke2] = ...
             set_diffusion_2a(grid, T_sflux1, T_sflux2, T_uflux1, T_uflux2, tke1, tke2, approx_d)

% Set diffusion coefficient at p-levels for diffusion
% of entropy, water, and w
% and at w-levels for diffusion of u, v, and tke

% This version uses turbulent KE and length scales to
% estimate diffusion coefficients

% Version 2a is equivalent to version 2 but is expressed in terms of input
% tke and time scales rather than tke and length scales; this make it
% easier to relate to MYNN coefficients and easier to bound rates / time
% scales to stabilize the MYNN scheme in growing turbulence

% The factor of 2/3 assumes M33 \approx q^2 / 3 = tke * 2/3

nz = grid.nz;

kdifft1 = T_sflux1.*tke1*2/3;
kdiffq1 = kdifft1;
kdiffw1 = T_uflux1.*tke1*2/3;
kdiffu1(2:nz) = grid.abovew(2:nz).*kdiffw1(2:nz) + grid.beloww(2:nz).*kdiffw1(1:nz-1);
kdiffu1(1)    = grid.extrapb1*kdiffw1(1) + grid.extrapb2*kdiffw1(2);
kdiffu1(nz+1) = grid.extraptnz*kdiffw1(nz) + grid.extraptnzm*kdiffw1(nz-1);
kdifftke1 = kdiffu1;


kdifft2 = T_sflux2.*tke2*2/3;
kdiffq2 = kdifft2;
kdiffw2 = T_uflux2.*tke2*2/3;
kdiffu2(2:nz) = grid.abovew(2:nz).*kdiffw2(2:nz) + grid.beloww(2:nz).*kdiffw2(1:nz-1);
kdiffu2(1)    = grid.extrapb1*kdiffw2(1) + grid.extrapb2*kdiffw2(2);
kdiffu2(nz+1) = grid.extraptnz*kdiffw2(nz) + grid.extraptnzm*kdiffw2(nz-1);
kdifftke2 = kdiffu2;

% ------

% Approximation d: no diffusive fluxes in fluid 2 (except surface flux)
if approx_d
    kdifft2 = 0*kdifft1;
    kdiffq2 = 0*kdiffq1;
    kdiffw2 = 0*kdiffw1;
    kdiffu2 = 0*kdiffu1;
    kdifftke2 = 0*kdifftke1;
    kdifft2(1) = kdifft1(1);
    kdiffq2(1) = kdiffq1(1);
    kdiffw2(1) = kdiffw1(1);
end


end