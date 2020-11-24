% Check that the increments of specific quantities are
% equal when the properties of the two fluids are equal, even
% if sigm1 and sigma2 depend on z.

% Call equate_fluid_properties before calling tendencies,
% then call this routine from quasinewton

% See also check_equal_tendencies

kkk = 1:4;

% Entropy
disp(' ')
disp('diff inc_eta')
inc_eta2(kkk) - inc_eta1(kkk)

% Water
disp(' ')
disp('diff inc_q')
inc_q2(kkk) - inc_q1(kkk)

% Vertical velocity
disp(' ')
disp('diff inc_w')
inc_w2(kkk) - inc_w1(kkk)

% Horizontal velocity
disp(' ')
disp('diff inc_u')
inc_u2(kkk) - inc_u1(kkk)
disp(' ')
disp('diff inc_v')
inc_v2(kkk) - inc_v1(kkk)

% TKE
disp(' ')
disp('diff inc_tke')
inc_tke2(kkk) - inc_tke1(kkk)

% eta variance
disp(' ')
disp('diff inc_vareta')
inc_vareta2(kkk) - inc_vareta1(kkk)

% q variance (advective form)
disp(' ')
disp('diff inc_varq')
inc_varq2(kkk) - inc_varq1(kkk)

pause
