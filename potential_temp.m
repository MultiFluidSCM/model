function [ theta ] = potential_temp( T, p, phys, therm )
% POTENTIAL_TEMP compute potential temperature given T and p

theta = T * (phys.p00/p)^(therm.kappa);

end

