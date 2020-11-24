function [ theta_rho ] = density_potential_temp( rho, p, phys, therm )
% DENSITY_POTENTIAL_TEMP compute density potential temperature given rho and p

theta_rho = phys.p00*(p/phys.p00)^(1 - therm.kappa) / (therm.Rd * rho);

end

