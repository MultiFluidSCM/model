function eos = find_eos_a(grid, state, constants)

% Find residuals in equations of state

% Also compute terms needed in linearization. The quantity needed
% is the derivative of (rho' / rho^2) wrt p', eta', and q'

% This version uses the mean fluid and fluid 2 thermodynamic equations
% and computes the implied rho^w and residuals for fluid 1

nz = grid.nz;
nzp = nz + 1;

% Mean fluid density, temperature, and moisture
m1  = state.fluid(1).m;
m2  = state.fluid(2).m;
rho = m1 + m2;
m1bar  = weight_to_w(grid,m1);
m2bar  = weight_to_w(grid,m2);
rhobar = weight_to_w(grid,rho);
T   = (m1.*state.fluid(1).T + m2.*state.fluid(2).T)./rho;
q   = (m1bar.*state.fluid(1).q   + m2bar.*state.fluid(2).q  )./rhobar;
eta = (m1bar.*state.fluid(1).eta + m2bar.*state.fluid(2).eta)./rhobar;
Tw  = (m1bar.*state.fluid(1).Tw  + m2bar.*state.fluid(2).Tw )./rhobar;

% On p levels
for k = 1:nz
    
    % Fluid 1
    qbar   = grid.abovep(k)*state.fluid(1).q(k+1) ...
           + grid.belowp(k)*state.fluid(1).q(k);
    etabar = grid.abovep(k)*state.fluid(1).eta(k+1) ...
           + grid.belowp(k)*state.fluid(1).eta(k);
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(state.p(k),          ...
                                                   state.fluid(1).T(k), ...
                                                   qbar,                ...
                                                   constants.therm);
    eos.rho1(k) = 1.0/gp;
    eos.sigma1(k) = state.fluid(1).m(k)*gp;
    % eos.res_sigma(k) = 1 - state.fluid(1).m(k)*gp;
    eos.drdp1(k)   = (gpt*gpt/gtt - gpp);
    eos.drdetap1(k) = gpt/gtt;
    eos.drdqp1(k)   = (gpt*gtw/gtt - gpw);
    
    % Fluid 2
    qbar   = grid.abovep(k)*state.fluid(2).q(k+1) ...
           + grid.belowp(k)*state.fluid(2).q(k);
    etabar = grid.abovep(k)*state.fluid(2).eta(k+1) ...
           + grid.belowp(k)*state.fluid(2).eta(k);
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(state.p(k),          ...
                                                   state.fluid(2).T(k), ...
                                                   qbar,                ...
                                                   constants.therm);
    eos.rho2(k) = 1.0/gp;
    eos.res_rho2(k) = 0.0;
    eos.sigma2(k) = state.fluid(2).m(k)*gp;
    % eos.res_sigma(k) = eos.res_sigma(k) - state.fluid(2).m(k)*gp;
    eos.res_etap2(k) = etabar + gt;
    eos.drdp2(k)   = (gpt*gpt/gtt - gpp);
    eos.drdetap2(k) = gpt/gtt;
    eos.drdqp2(k)   = (gpt*gtw/gtt - gpw);
    
    % Mean fluid
    qbar   = grid.abovep(k)*q(k+1) ...
           + grid.belowp(k)*q(k);
    etabar = grid.abovep(k)*eta(k+1) ...
           + grid.belowp(k)*eta(k);
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(state.p(k),          ...
                                                   T(k),                ...
                                                   qbar,                ...
                                                   constants.therm);
    res_rho = 1.0/rho(k) - gp;
    res_etap = etabar + gt;
    
    % Residual in sigma equation
    eos.res_sigma(k) = 1 - eos.sigma1(k) - eos.sigma2(k);

    % Implied fluid 1 fields
    eos.res_rho1(k)  = (rho(k)*res_rho  - m2(k)*eos.res_rho2(k) - eos.res_sigma(k))/m1(k);
    eos.res_etap1(k) = (rho(k)*res_etap - m2(k)*eos.res_etap2(k))/m1(k);
    
end


% On w levels
for k = 1:nzp
    
    if k == 1
        pbar = grid.extrapb1*state.p(1) + grid.extrapb2*state.p(2);
    elseif k == nzp
        pbar = grid.extraptnz*state.p(nz) + grid.extraptnzm*state.p(nz-1);
    else
        pbar   = grid.abovew(k)*state.p(k) ...
               + grid.beloww(k)*state.p(k-1);
    end
    
    % Fluid2
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(pbar,                 ...
                                                   state.fluid(2).Tw(k), ...
                                                   state.fluid(2).q(k),  ...
                                                   constants.therm);
    eos.rhow2(k) = 1.0/gp;
    eos.res_eta2(k) = state.fluid(2).eta(k) + gt;
    eos.drdpbar2(k) = (gpt*gpt/gtt - gpp);
    eos.drdeta2(k)   = gpt/gtt;
    eos.drdq2(k)    = (gpt*gtw/gtt - gpw);
    
     % Mean fluid
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(pbar,                 ...
                                                   Tw(k),                ...
                                                   q(k),                 ...
                                                   constants.therm);
    rhow = 1.0/gp;
    res_eta = eta(k) + gt;
    
    % Implied fluid 1 fields
    eos.rhow1(k) = m1bar(k)/(rhobar(k)/rhow - m2bar(k)/eos.rhow2(k));
    eos.res_eta1(k) = (rhobar(k)*res_eta - m2bar(k)*eos.res_eta2(k))/m1bar(k);
    eos.drdpbar1(k) = (gpt*gpt/gtt - gpp);
    eos.drdeta1(k)   = gpt/gtt;
    eos.drdq1(k)    = (gpt*gtw/gtt - gpw);
    
end

disp('** check this is correct for switch a **')
pause

% Vertical pressure gradient
dpdz(2:nz)   = (state.p(2:nz) - state.p(1:nz-1))./grid.dzw(2:nz);
dpdz(1) = dpdz(2);
dpdz(nzp) = dpdz(nz);
dpdzbar = grid.abovep.*dpdz(2:nzp) + grid.belowp.*dpdz(1:nz);

% Buoyancy frequency squared
% Note there are several ways to define this, depending on what is assumed
% about the parcel and the environment through which it moves. Here we
% consider  a parcel that conserves its eta and q moving relative to a
% background in which eta1, eta2, q1, q2, and p are steady.

deta1dz = (state.fluid(1).eta(2:nzp) - state.fluid(1).eta(1:nz))./grid.dzp;
dq1dz = (state.fluid(1).q(2:nzp) - state.fluid(1).q(1:nz))./grid.dzp;
deta2dz = (state.fluid(2).eta(2:nzp) - state.fluid(2).eta(1:nz))./grid.dzp;
dq2dz = (state.fluid(2).q(2:nzp) - state.fluid(2).q(1:nz))./grid.dzp;
env_term = eos.sigma1.*(eos.drdetap1.*deta1dz + eos.drdqp1.*dq1dz + eos.drdp1.*dpdzbar) ...
         + eos.sigma2.*(eos.drdetap2.*deta2dz + eos.drdqp2.*dq2dz + eos.drdp2.*dpdzbar);
eos.nsq1 = -constants.phys.gravity*eos.rho1.*(env_term - eos.drdp1.*dpdzbar);
eos.nsq2 = -constants.phys.gravity*eos.rho2.*(env_term - eos.drdp2.*dpdzbar);



end

