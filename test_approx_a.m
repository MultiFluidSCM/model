% Test implementation of approximation a by checking how well the equation
% of state is satisfied
% 
% Approximation a: Use mean fluid and fluid 2 thermodynamic equations
eos = find_eos_a(grid, state_new, constants);
%
disp('Appriximation (a)')
disp(['Res rho1  = ', num2str(eos.res_rho1(1:6))])
disp(['Res etap1 = ', num2str(eos.res_etap1(1:6))])
disp(['Res eta1  = ', num2str(eos.res_eta1(1:6))])

% Use fluid 1 and fluid 2 thermodynamic equations
eos = find_eos(grid, state_new, constants);
disp('No appriximation (a)')
disp(['Res rho1  = ', num2str(eos.res_rho1(1:6))])
disp(['Res etap1 = ', num2str(eos.res_etap1(1:6))])
disp(['Res eta1  = ', num2str(eos.res_eta1(1:6))])