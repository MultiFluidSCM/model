% Equate all the properties of the two fluids except the sigmas

m_temp = state_old.fluid(2).m;
state_old.fluid(2) = state_old.fluid(1);
state_old.fluid(1).w(:) = 0;
state_old.fluid(2).w(:) = 0;
state_old.fluid(2).m = m_temp;

% rho = state_old.fluid(1).m + state_old.fluid(2).m;
% state_old.fluid(1).m = 0.5*rho;
% state_old.fluid(2).m = 0.5*rho;