function [kdifft1,kdiffq1,kdiffw1,kdifft2,kdiffq2,kdiffw2,kdiffu1,kdiffu2,kdifftke1,kdifftke2] = ...
             set_diffusion(zp, zw, zstar, wstar, ustar_by_wstar, k_vonkarman, zrough, approx_d)

% Set diffusion coefficient at p-levels for diffusion
% of entropy, water, and w

% First calculate dimensionless height at p-levels
zhat = min(1,(zp + zrough)/(1.1*zstar + zrough));

kscale = max(zstar*wstar,1);

% Now the diffusion coefficients
kdifft1 = kscale*k_vonkarman...
         *(((ustar_by_wstar)^3 + 39*k_vonkarman*zhat).^(1/3))...
         .*zhat...
         .*((1 - zhat).^2);
kdiffq1 = kdifft1;
kdiffw1 = kdifft1;


% Set diffusion coefficient at w-levels for diffusion
% of u and v and tke

% First calculate dimensionless height at w-levels
zhat = min(1,(zw + zrough)/(1.1*zstar + zrough));

% Now the diffusion coefficients
kdiffu1 = kscale*k_vonkarman...
         *(((ustar_by_wstar)^3 + 39*k_vonkarman*zhat).^(1/3))...
         .*zhat...
         .*((1 - zhat).^2);
kdifftke1 = kdiffu1;


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
else
    kdifft2 = kdifft1;
    kdiffq2 = kdiffq1;
    kdiffw2 = kdiffw1;
    kdiffu2 = kdiffu1;
    kdifftke2 = kdifftke1;
end


end