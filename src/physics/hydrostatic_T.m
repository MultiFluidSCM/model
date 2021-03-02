% Set up profiles in hydrostatic balance with a specified T and q

% Surface temperature and water
tt = initial_T(zw(1));
ww = initial_q(zw(1));

% Other surface fields
[g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(p0s,tt,ww,constants.therm);

% Save surface eta, water, and T
eta(1) = -gt;
water(1) = ww;
Tw(1) = tt;


% Given T at zw(2) ...
tt = initial_T(zw(2));
ww = initial_q(zw(2));

% ... determine pressure at level 1, such that when extrapolated to the surface
% we get the desired surface pressure.
% First guess
p1 = p0s;
p2 = p1;
for iter1 = 1:10

  % Find p2 that gives hydrostatic balance
  for iter2 = 1:10
    pbar = abovew(2)*p2 + beloww(2)*p1;
    dpbydz = (p2 - p1)/dzw(2);
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(pbar,tt,ww,constants.therm);
    r1 = gp*dpbydz + gravity;
    % Newton update
    aa = abovew(2)*gpp*dpbydz + gp/dzw(2);
    p_inc = -r1/aa;
    p2 = p2 + p_inc;
  end

  % Estimated surface pressure
  psurf = grid.extrapb2*p2 + grid.extrapb1*p1;
  
  % New value of p  
  p1 = p1 + p0s - psurf;

end

% Save p
p(1) = p1;


% Integrate the hydrostatic relation to find the pressure
for k = 2:nz
  p2 = p1;
  % Given T at zw(k) ...
  tt = initial_T(zw(k));
  ww = initial_q(zw(k));
  % ... find p2 and T that give hydrostatic balance with the right eta
  for iter2 = 1:10
    pbar = abovew(k)*p2 + beloww(k)*p1;
    dpbydz = (p2 - p1)/dzw(k);
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(pbar,tt,ww,constants.therm);
    r1 = gp*dpbydz + gravity;
    % Newton update
    aa = abovew(k)*gpp*dpbydz + gp/dzw(k);
    p_inc = -r1/aa;
    p2 = p2 + p_inc;
  end
  p1 = p2;
  
  % Save eta-level values
  eta(k) = -gt;
  water(k) = ww;
  Tw(k) = tt;
  
  % And rho-level p
  p(k) = p2;
  
end

% Top boundary values
tt = initial_T(zw(nzp));
ww = initial_q(zw(nzp));
pbar = grid.extraptnz*p(nz) + grid.extraptnzm*p(nz-1);
[g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(pbar,tt,ww,constants.therm);
eta(nzp) = -gt;
water(nzp) = ww;
Tw(nzp) = tt;

