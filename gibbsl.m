function [ g, gp, gt, gpp, gpt, gtt ] = gibbsl( p, T, therm )
%GIBBSL Gibbs function and derivatives for liquid water
%   Given pressure p and temperature T, compute the Gibbs
%   thermodynamic potential for liquid water, along with its first
%   and second derivatives.

%   This version assumes a perfect incompressible liquid.

%   The latent heating term has been moved from vapour to liquid cf T17

% Unpack constants for clarity
T0     = therm.T0;
Cl     = therm.Cl;
alpha0 = therm.alpha0;
p0v    = therm.p0v;
L00    = therm.L00;

% Useful intermediate quantities
LT = log(T/T0);

% Gibbs function
g = -Cl*LT*T + alpha0*(p - p0v*T/T0) - L00*(1 - T/T0);

% First derivatives
gp = alpha0;
gt = -Cl*(1 + LT) - alpha0*p0v/T0 + L00/T0;

% Second derivatives
gpp = 0;
gpt = 0;
gtt = -Cl/T;


end

