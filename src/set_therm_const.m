function therm = set_therm_const( )

% Set thermodynamic constants

therm.Rd = 287.00;         % Gas constant for dry air
therm.Rv = 461.00;         % Gas constant for water vapour
therm.epsilon = therm.Rd/therm.Rv;
therm.Cpd = 1004.0;        % Specific heat capacity at constant pressure, dry air
therm.Cpv = 1885.0;        % Specific heat capacity at constant pressure, vapour
therm.Cl = 4186.0;         % Specific heat capacity of liquid water
therm.Cf = 2106.0;         % Specific heat capacity of frozen water
therm.Cvd = therm.Cpd - therm.Rd;
therm.Cvv = therm.Cpv - therm.Rv;
therm.kappa = therm.Rd/therm.Cpd;
therm.L0 = 2.501e6;        % Specific latent heat of vaporization at T0 and p0d [p0v]
therm.L0s = 2.834e6;       % Specific latent heat of sublimation at T0 and p0d [p0v]
                           % L0f = L0s - L0 = 0.333e6
therm.alpha0 = 0.001;      % Specific volume of liquid water
therm.alpha0f = 0.0011;    % Specific volume of frozen water

therm.p0d = 100000.0;      % Reference pressure for dry air
therm.p0v = 611.657;       % Reference pressure for vapour, equal to saturation
                           % vapour pressure at T0
therm.p0vf = 611.657;      % Reference pressure for vapour over ice, equal to saturation
                           % vapour pressure at T0
therm.T0 = 273.16;         % Reference temperature
therm.L00 = therm.L0 - (therm.Cpv - therm.Cl)*therm.T0 + therm.alpha0*therm.p0v; % Specific latent heat of vaporization
                           % extrapolated to 0K and p = 0
therm.L00s = therm.L0s - (therm.Cpv - therm.Cf)*therm.T0 + therm.alpha0f*therm.p0v; % Specific latent heat of sublimation
                           % extrapolated to 0K and p = 0

end

