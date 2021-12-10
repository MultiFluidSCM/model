function [ KHdrag ] = find_KHdrag( m1, m2, m1bar, m2bar, w1, w2, tke1, tke2, scales, grid )

%
% Estimate effective drag due to two-fluid Kelvin-Helmholtz instability.
% The drag will provide a sink of resolved KE and hence a source of tke.
% It is computed such that tke equilibrates when eddy diffusivity is strong
% enough to controll the instability.

% Here we just consider the vertical component

% Unpack fields for clarity of code
LT1 = scales.L_turb1;
LT2 = scales.L_turb2;


% Which direction is the drag?
s = sign(w2 - w1);

% Assumed length scale of instability
% *** will need to rethink this in the cloud layer ***
LKH = scales.zstar;

% Sum of tke divided by LT^2
% ebyL2 = m1.*tke1./(LT1.*LT1) + m2.*tke2./(LT2.*LT2);
ebyL2 = (m1.*tke1 + m2.*tke2)/(LKH*LKH);

% Interpolate to w-levels
ebyL2w = weight_to_w(grid,ebyL2);

% Calculate drag
KHdrag = LKH * s .* (sqrt(m1bar.*m2bar)./(m1bar + m2bar)) .* ebyL2w;


end

