function [ g, gp, gt, gpp, gpt, gtt ] = gibbsf( p, T )
%GIBBSF Gibbs function and derivatives for frozen water
%   Given pressure p and temperature T, compute the Gibbs
%   thermodynamic potential for frozen water, along with its first
%   and second derivatives.

%   This version assumes a perfect incompressible solid.

%   The latent heating term has been moved from vapour to liquid cf T17

% Set constants
constants

% Useful intermediate quantities
LT = log(T/T0);

% Gibbs function
g = -Cf*LT*T + alpha0f*(p - p0vf*T/T0) ...
   - L00s*(1 - T/T0);

% First derivatives
gp = alpha0f;
gt = -Cf*(1 + LT) - alpha0f*p0vf/T0 ...
    + L00s/T0;

% Second derivatives
gpp = 0;
gpt = 0;
gtt = -Cf/T;


end

