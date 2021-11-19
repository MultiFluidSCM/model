% Set up profiles in hydrostatic balance with a specified thetal and q

% Determine surface entropy
thetal = initial_theta(zw(1), settings.initial_theta);
ww = initial_q(zw(1), settings.initial_qv);
eta00 = thetal2eta(thetal,ww,constants.therm,constants.phys.p00);

% Surface exner (for estimating surface temperature when p0s differs from p00)
exner0s = (p0s/constants.phys.p00)^constants.therm.kappa;





% Search by bisection to get a sensible starting value of
% surface temperature
t2 = (thetal + 1.0)*exner0s;
t1 = (t2 - 4.0d3*ww - 1.0)*exner0s;
[g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(p0s,t1,ww,constants.therm);
r1 = gt + eta00;
[g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(p0s,t2,ww,constants.therm);
r2 = gt + eta00;

if (r1*r2 > 0.0d0)
  disp('Need wider bounds to start bisection search')
  pause
end
for iter1 = 1:5
  tt = 0.5d0*(t1 + t2);
  [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(p0s,tt,ww,constants.therm);
  rr = gt + eta00;
  if (r1*rr > 0.0d0)
    t1 = tt;
    r1 = rr;
  else
    t2 = tt;
    r2 = rr;
  end
end

% And finish off with Newton iteration
for iter1 = 1:5
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(p0s,tt,ww,constants.therm);
    r1 = gt + eta00;
    t_inc = - r1/gtt;
    tt = tt + t_inc;
end

% Save surface eta, water, and T
eta(1) = eta00;
water(1) = ww;
Tw(1) = tt;

% Given theta at zw(2) ...
ww = initial_q(zw(2), settings.initial_qv);
thetal = initial_theta(zw(2), settings.initial_theta);
eta00 = thetal2eta(thetal,ww,constants.therm,constants.phys.p00);

% ... determine pressure at level 1, such that when extrapolated to the surface
% we get the desired surface pressure.
% First guess
p1 = p0s;
p2 = p1;
for iter1 = 1:10

  % Find p2 and T that give hydrostatic balance with the right eta
  for iter2 = 1:10
    pbar = abovew(2)*p2 + beloww(2)*p1;
    dpbydz = (p2 - p1)/dzw(2);
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(pbar,tt,ww,constants.therm);
    r1 = gp*dpbydz + gravity;
    r2 = gt + eta00;
    % Newton update
    aa = abovew(2)*gpp*dpbydz + gp/dzw(2);
    bb = gpt;
    cc = abovew(2)*gpt;
    dd = gtt;
    rdet = 1.0d0/(aa*dd - bb*cc);
    p_inc = rdet*(bb*r2 - dd*r1);
    t_inc = rdet*(cc*r1 - aa*r2);
    p2 = p2 + p_inc;
    tt = tt + t_inc;
  end

  % Estimated surface pressure
  psurf = grid.extrapb2*p2 + grid.extrapb1*p1;
  
  % New value of p  
  p1 = p1 + p0s - psurf;

end

% Save rho-level pressure
p(1) = p1;


% Integrate the hydrostatic relation to find the pressure
for k = 2:nz
  p2 = p1;
  
  % Given theta at zw(k) ...
  ww = initial_q(zw(k), settings.initial_qv);
  thetal = initial_theta(zw(k), settings.initial_theta);
  eta00 = thetal2eta(thetal,ww,constants.therm,constants.phys.p00);
  % ... find p2 and T that give hydrostatic balance with the right eta
  for iter2 = 1:10
    pbar = abovew(k)*p2 + beloww(k)*p1;
    dpbydz = (p2 - p1)/dzw(k);
    
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(pbar,tt,ww,constants.therm);
    r1 = gp*dpbydz + gravity;
    r2 = gt + eta00;
    % Newton update
    aa = 0.5d0*gpp*dpbydz + gp/dzw(k);
    bb = gpt;
    cc = 0.5d0*gpt;
    dd = gtt;
    rdet = 1.0d0/(aa*dd - bb*cc);
    p_inc = rdet*(bb*r2 - dd*r1);
    t_inc = rdet*(cc*r1 - aa*r2);
    p2 = p2 + p_inc;
    tt = tt + t_inc;
  end
  p1 = p2;
  
  % Save eta-level values
  eta(k) = eta00;
  water(k) = ww;
  Tw(k) = tt;
  
  % And rho-level p
  p(k) = p2;
  
end

% Top boundary values
ww = initial_q(zw(nzp), settings.initial_qv);
thetal = initial_theta(zw(nzp), settings.initial_theta);
eta00 = thetal2eta(thetal,ww,constants.therm,constants.phys.p00);
eta(nzp) = eta00;
water(nzp) = ww;
pbar = grid.extraptnz*p(nz) + grid.extraptnzm*p(nz-1);
tt = Tw(nz);
for iter = 1:5
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(pbar,tt,ww,constants.therm);
    r1 = gt + eta00;
    t_inc = -r1/gtt;
    tt = tt + t_inc;
end
Tw(nzp) = tt;

