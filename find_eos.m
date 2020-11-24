function eos = find_eos(grid, state, constants)

% Determine densities and find residuals in equations of state

% Also compute terms needed in linearization. The quantity needed
% is the derivative of (rho' / rho^2) wrt p', eta', and q'

nz = grid.nz;
nzp = nz + 1;

% On p levels
for k = 1:nz
    
    % Fluid 1
    qbar   = grid.aboves(k)*state.fluid(1).q(k+1) ...
           + grid.belows(k)*state.fluid(1).q(k);
    etabar = grid.aboves(k)*state.fluid(1).eta(k+1) ...
           + grid.belows(k)*state.fluid(1).eta(k);
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(state.p(k),          ...
                                                   state.fluid(1).T(k), ...
                                                   qbar,                ...
                                                   constants.therm);
    eos.rho1(k) = 1.0/gp;
    eos.res_rho1(k) = 0.0;
    eos.sigma1(k) = state.fluid(1).m(k)*gp;
    % eos.res_sigma(k) = 1 - state.fluid(1).m(k)*gp;
    eos.res_etap1(k) = etabar + gt;
    eos.drdp1(k)    = (gpt*gpt/gtt - gpp);
    eos.drdetap1(k) =  gpt/gtt;
    eos.drdqp1(k)   = (gpt*gtw/gtt - gpw);
    
    % Fluid 2
    qbar   = grid.aboves(k)*state.fluid(2).q(k+1) ...
           + grid.belows(k)*state.fluid(2).q(k);
    etabar = grid.aboves(k)*state.fluid(2).eta(k+1) ...
           + grid.belows(k)*state.fluid(2).eta(k);
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(state.p(k),          ...
                                                   state.fluid(2).T(k), ...
                                                   qbar,                ...
                                                   constants.therm);
    eos.rho2(k) = 1.0/gp;
    eos.res_rho2(k) = 0.0;
    eos.sigma2(k) = state.fluid(2).m(k)*gp;
    % eos.res_sigma(k) = eos.res_sigma(k) - state.fluid(2).m(k)*gp;
    eos.res_etap2(k) = etabar + gt;
    eos.drdp2(k)    = (gpt*gpt/gtt - gpp);
    eos.drdetap2(k) =  gpt/gtt;
    eos.drdqp2(k)   = (gpt*gtw/gtt - gpw);

    % Residual in sigma equation
    eos.res_sigma(k) = 1 - eos.sigma1(k) - eos.sigma2(k);

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
    
    % Fluid 1
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(pbar,                 ...
                                                   state.fluid(1).Tw(k), ...
                                                   state.fluid(1).q(k),  ...
                                                   constants.therm);
    eos.rhow1(k) = 1.0/gp;
    eos.res_eta1(k) = state.fluid(1).eta(k) + gt;
    eos.drdpbar1(k) = (gpt*gpt/gtt - gpp);
    eos.drdeta1(k)  =  gpt/gtt;
    eos.drdq1(k)    = (gpt*gtw/gtt - gpw);
  
    % Fluid 2
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(pbar,                 ...
                                                   state.fluid(2).Tw(k), ...
                                                   state.fluid(2).q(k),  ...
                                                   constants.therm);
    eos.rhow2(k) = 1.0/gp;
    eos.res_eta2(k) = state.fluid(2).eta(k) + gt;
    eos.drdpbar2(k) = (gpt*gpt/gtt - gpp);
    eos.drdeta2(k)  =  gpt/gtt;
    eos.drdq2(k)    = (gpt*gtw/gtt - gpw);

    eos.pbar(k) = pbar;
    
    % eos.theta1(k) = potential_temp(state.fluid(1).Tw(k),pbar,constants.phys,constants.therm);
    % eos.theta2(k) = potential_temp(state.fluid(2).Tw(k),pbar,constants.phys,constants.therm);
    eos.theta1(k) = eta2thetal( state.fluid(1).eta(k), ...
                                state.fluid(1).q(k), ...
                                state.fluid(1).Tw(k), ...
                                constants.therm, constants.phys.p00);
    eos.theta2(k) = eta2thetal( state.fluid(2).eta(k), ...
                                state.fluid(2).q(k), ...
                                state.fluid(2).Tw(k), ...
                                constants.therm, constants.phys.p00);
    eos.theta_rho1(k) = density_potential_temp(eos.rhow1(k),pbar,constants.phys,constants.therm);
    eos.theta_rho2(k) = density_potential_temp(eos.rhow2(k),pbar,constants.phys,constants.therm);

    eos.rho_ptl1(k) = potential_density(constants.phys.p00,state.fluid(1).eta(k),...
                                        state.fluid(1).q(k),state.fluid(1).Tw(k),constants.therm);
    eos.rho_ptl2(k) = potential_density(constants.phys.p00,state.fluid(2).eta(k),...
                                        state.fluid(2).q(k),state.fluid(2).Tw(k),constants.therm);
    
end


% Vertical pressure gradient
dpdz(2:nz)   = (state.p(2:nz) - state.p(1:nz-1))./grid.dzw(2:nz);
dpdz(1) = dpdz(2);
dpdz(nzp) = dpdz(nz);
dpdzbar = grid.abovep.*dpdz(2:nzp) + grid.belowp.*dpdz(1:nz);

% Buoyancy frequency squared
% Note there are several ways to define this, depending on what is assumed
% about the parcel and the environment through which it moves. Here we
% consider a parcel that conserves its eta and q moving relative to a
% background in which eta1, eta2, q1, q2, and p are steady.

% Old calculation of environment term - very sensitive to strong gradients
% in the lowest layer
%deta1dz = (state.fluid(1).eta(2:nzp) - state.fluid(1).eta(1:nz))./grid.dzp;
%dq1dz = (state.fluid(1).q(2:nzp) - state.fluid(1).q(1:nz))./grid.dzp;
%deta2dz = (state.fluid(2).eta(2:nzp) - state.fluid(2).eta(1:nz))./grid.dzp;
%dq2dz = (state.fluid(2).q(2:nzp) - state.fluid(2).q(1:nz))./grid.dzp;
%env_term = eos.sigma1.*(eos.drdetap1.*deta1dz + eos.drdqp1.*dq1dz + eos.drdp1.*dpdzbar) ...
%         + eos.sigma2.*(eos.drdetap2.*deta2dz + eos.drdqp2.*dq2dz + eos.drdp2.*dpdzbar);
%eos.nsq1x = -constants.phys.gravity*eos.rho1.*(env_term - eos.drdp1.*dpdzbar);
%eos.nsq2x = -constants.phys.gravity*eos.rho2.*(env_term - eos.drdp2.*dpdzbar);

drho1dz = (eos.rhow1(2:nzp) - eos.rhow1(1:nz))./grid.dzp;
drho2dz = (eos.rhow2(2:nzp) - eos.rhow2(1:nz))./grid.dzp;
drhodz = eos.sigma1.*drho1dz + eos.sigma2.*drho2dz;

rho = state.fluid(1).m + state.fluid(2).m;
env_term_alt = drhodz./(rho.*rho);
eos.nsq1 = -constants.phys.gravity*eos.rho1.*(env_term_alt - eos.drdp1.*dpdzbar);
eos.nsq2 = -constants.phys.gravity*eos.rho2.*(env_term_alt - eos.drdp2.*dpdzbar);




end