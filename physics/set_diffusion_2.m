function [kdifft1,kdiffq1,kdiffw1,kdifft2,kdiffq2,kdiffw2,kdiffu1,kdiffu2,kdifftke1,kdifftke2] = ...
             set_diffusion_2(grid, L_turb1, L_turb2, tke1, tke2, approx_d)

% Set diffusion coefficient at p-levels for diffusion
% of entropy, water, and w
% and at w-levels for diffusion of u, v, and tke

% This version uses turbulent KE and length scales to
% estimate diffusion coefficients


nz = grid.nz;

kdifft1 = L_turb1.*sqrt(tke1);
kdiffq1 = kdifft1;
kdiffw1 = kdifft1;
kdiffu1(2:nz) = grid.abovew(2:nz).*kdifft1(2:nz) + grid.beloww(2:nz).*kdifft1(1:nz-1);
kdiffu1(1)    = grid.extrapb1*kdifft1(1) + grid.extrapb2*kdifft1(2);
kdiffu1(nz+1) = grid.extraptnz*kdifft1(nz) + grid.extraptnzm*kdifft1(nz-1);
kdifftke1 = kdiffu1;


kdifft2 = L_turb2.*sqrt(tke2);
kdiffq2 = kdifft2;
kdiffw2 = kdifft2;
kdiffu2(2:nz) = grid.abovew(2:nz).*kdifft2(2:nz) + grid.beloww(2:nz).*kdifft2(1:nz-1);
kdiffu2(1)    = grid.extrapb1*kdifft2(1) + grid.extrapb2*kdifft2(2);
kdiffu2(nz+1) = grid.extraptnz*kdifft2(nz) + grid.extraptnzm*kdifft2(nz-1);
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