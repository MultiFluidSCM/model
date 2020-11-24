function [ rho_ptl ] = potential_density( p00, eta, q, T, therm )
%POTENTIAL_DENSITY Compute potential density
%   Compute the density that this parcel would have if taken to pressure
%   p00. The parcel temperature T provides a first guess for an iterative
%   solution.

for iter = 1:3
    
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(p00,T,q,therm);
    res = gt + eta;
    T_inc = -res/gtt;
    T = T + T_inc;

end

rho_ptl = 1.0/gp;
%pause