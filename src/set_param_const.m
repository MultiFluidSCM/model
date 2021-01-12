function param = set_param_const( )

% Set parameterization constants

% Note: setting bentrain < 1 introduces a positive feedback that can
% cause the solver to diverge when the entrainment timescale is short
% (e.g. idealised dry CBL case with large initial surface heat flux).
% It appears to be less of a problem in more realistic cases that spin
% up more gradually.

param.sigma00   = 0.01;   % Background sigma2 when nothing is going on
% param.centrain  = 0.4;     % Entrainment coefficient
% param.cdetrain  = 0.7;     % Detrainment coefficient
param.bentrainw = 0.5;     % Factor for entrainment of w
param.bentraint = 1.0;     % Factor for entrainment of eta
param.bentrainq = 1.0;     % Factor for entrainment of water
param.bentrainu = 1.0;     % Factor for detrainment of u and v
param.bdetrainw = 1.0;     % Factor for detrainment of w
param.bdetraint = 1.0;     % Factor for detrainment of eta
param.bdetrainq = 1.0;     % Factor for detrainment of water
param.bdetrainu = 1.0;     % Factor for detrainment of u and v
param.confrac   = 0.10;    % Reference updraft mass fraction
%param.alpha_plume = 1.5;   % Constant for updraft eta and q contrast
param.zrough    = 0.1;     % Roughness length

% Magic numbers - dimensional constants that are not constants
% of nature - to be deprecated and avoided if at all possible
param.tke_min   = 1e-4 ;   % Minimum permitted tke (J / kg)
param.zstar_min = 50;      % Minimum allowed boundary layer depth (m)

end

