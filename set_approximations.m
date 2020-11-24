function [ switches ] = set_approximations( )
%SET_APPROXIMATIONS Decide which approximations to apply
% They are all controlled by switches

% Switch (a)
% 0 use fluid 1 thermodynamic equation;
% 1 replace fluid 1 thermodynamic equation by mean fluid
% thermodynamic equation
switches.a = 0;

% Switch (b)
% 0 use sum of fluid 1 and fluid 2 subfilter-scale fluxes in the
%   implied mean fluid equations
% 1 use a mean fluid subfilter-scale flux in the implied mean
%   fluid equations
switches.b = 0;

% Switch (c)
% 0 use full prognostic equations for fluid 2
% 1 neglect time derivatives in fluid 2 but retain full time
%   derivatives in mean fluid equations
switches.c = 0;

% Switch (d)
% 0 retain subfilter-scale fluxes in fluid 2
% 1 neglect subfilter-scale fluxes in fluid 2 except at first p level
switches.d = 0;

% Switch (e)
% 0 retain the full conditional filtering fluxes in the implied mean
%   fluid equation
% 1 make the small area approximation in the conditional filtering
%   fluxes in the implied mean fluid equation
switches.e = 0;

end

