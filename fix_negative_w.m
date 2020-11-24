function [inc_w1,inc_w2,inc_eta1,inc_eta2,inc_q1,inc_q2] = fix_negative_w( grid, state_in )
%FIXER Fix descent in updraft
%   Fix descent in updraft by homoegizing w, eta and q
%   between fluids

%  This version computes advective form increments allowing budgets
%  to be computed

% Remap mass to w levels
m1bar = weight_to_w(grid,state_in.fluid(1).m);
m2bar = weight_to_w(grid,state_in.fluid(2).m);

% Total mass at w levels
mtotbar = m1bar + m2bar;

% Is velocity difference of the wrong sign
lfix = 0.5*(1 + sign(state_in.fluid(1).w - state_in.fluid(2).w));
lfix(1) = 0;
lfix(end) = 0;

% Homogenize at levels where needed
w_homog   = (m1bar.*state_in.fluid(1).w   + m2bar.*state_in.fluid(2).w  )./mtotbar;
eta_homog = (m1bar.*state_in.fluid(1).eta + m2bar.*state_in.fluid(2).eta)./mtotbar;
q_homog   = (m1bar.*state_in.fluid(1).q   + m2bar.*state_in.fluid(2).q  )./mtotbar;

% And compute increments
inc_w1   = lfix.*(w_homog   - state_in.fluid(1).w);
inc_w2   = lfix.*(w_homog   - state_in.fluid(2).w);
inc_eta1 = lfix.*(eta_homog - state_in.fluid(1).eta);
inc_eta2 = lfix.*(eta_homog - state_in.fluid(2).eta);
inc_q1   = lfix.*(q_homog   - state_in.fluid(1).q);
inc_q2   = lfix.*(q_homog   - state_in.fluid(2).q);


end

