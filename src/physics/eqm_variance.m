function [ var1, var2 ] = eqm_variance( grid, q1, q2, m1bar, m2bar, M12bar, M21bar, ...
                                        rate1, rate2, D1, D2, corrde, qhat12, qhat21 )
%EQM_VARIANCE Estimate the SF-scale variance of a scalar
%   Estimate the SF-scale variance of a scalar field q by
%   assuming local equilibrium between sources and sinks.
%   This is the multi-fluid extension of MYNN level 2.5
%   for scalars. rate1 and rate2 are turbulent dissipation rates.
%   As an approximation, it is assumed that
%   entrained and detrained values of variance are `upwind' values
%   Note that if qhat's depend on variance then an iterative solution
%   is required.


nz = grid.nz;
nzp = nz + 1;

% Vertical mean of total density
rhobar = m1bar + m2bar;

% Vertical derivatives
dq1dz = (q1(2:nzp) - q1(1:nz))./grid.dzp;
dq2dz = (q2(2:nzp) - q2(1:nz))./grid.dzp;

% Downgradient flux source terms ...
dg1 = D1.*dq1dz;
dg2 = D2.*dq2dz;

% ... mapped to w levels
dg1w = weight_to_w(grid,dg1);
dg2w = weight_to_w(grid,dg2);

% `diffent' source terms
temp = (q2 - q1).*corrde./rhobar;
de1 = m1bar.*temp;
de2 = m2bar.*temp;

% Entrainment and detrainment source terms
Sd1 = M12bar.*(  q1.*q1 + q2.*q2 - 2*q1.*qhat12);
Se1 = M21bar.*(         2*q2.*q2 - 2*q1.*qhat21);
Se2 = M21bar.*(  q1.*q1 + q2.*q2 - 2*q2.*qhat21);
Sd2 = M12bar.*(2*q1.*q1          - 2*q2.*qhat12);


% Loop over levels
for k = 2:nz
    
    % At each level solve a 2x2 system for the variances
    
    % Right hand sides
    R1 = -2*dg1(k) - 2*de1(k) + Sd1(k) - Se1(k);
    R2 = -2*dg2(k) - 2*de2(k) + Se2(k) - Sd1(k);
    
    % Coefficients in linear system
    A11 =  M21bar(k) + 2*m1bar(k)*rate1(k);
    A12 = -M12bar(k);
    A21 = -M21bar(k);
    A22 =  M12bar(k) + 2*m2bar(k)*rate2(k);
    
    % And solve linear system
    rdet = 1/(A11*A22 - A12*A21);
    var1(k) = rdet*(A22*R1 - A12*R2);
    var2(k) = rdet*(A11*R2 - A21*R1);
    
end

% Crudely extrapolate to boundaries
var1(1)   = var1(2);
var1(nzp) = var1(nz);
var2(1)   = var2(2);
var2(nzp) = var2(nz);

