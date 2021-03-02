function [ LNB ] = LNBgas( eta, q, T, zw, rho, p, therm )
%LNBGAS Find the level of neutral buoyancy of a surface parcel
%       ignoring condensation and latent heat release


% eta and q are properties of the parcel to be tested
% T is a first guess for the parcel temperature for the
% iterative solver
a = 1 - q;

% Initialize search
k = 1;
klast = numel(zw);
buoy_test = 0;
buoyant = 1;

% Search for the level at which the parcel is no longer buoyant
while buoyant & k < klast
    
    % Next level
    buoy_prev = buoy_test;
    k = k + 1;
    
    % Find the density of the lifted parcel at the pressure of level k
    % First need to solve for the temperature - Newton iteration
    % A single iteration might be enough, but let's use 2
    for iter = 1:2
        % Ignore condensation so use gibbs for air plus vapour
        [g,gp,gt,ga,gpp,gpt,gpa,gtt,gta,gaa,gtwv] = gibbsav(p(k),T,a,therm);
        res = gt + eta;
        T_inc = -res/gtt;
        T = T + T_inc;
    end
    buoy_test = rho(k) - 1/gp;
    buoyant = buoy_test > 0;

end

% Interpolate between last two values to get a more accurate estimate
LNB = (zw(k-1)*buoy_test - zw(k)*buoy_prev)/(buoy_test - buoy_prev);


