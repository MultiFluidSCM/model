function [ tend, relabel, eos, force, scales, surface_flux, budgets, work ] ...
        = tendencies( grid, state, settings, t, dt, switches , old_diff)

% Compute time tendencies of all fields for a given state

% ------

% Unpack fields to improve code clarity
nz = grid.nz;
nzp = nz + 1;
dzp = grid.dzp;
dzw = grid.dzw;
abovep = grid.abovep;
belowp = grid.belowp;
p    = state.p;
m1   = state.fluid(1).m;
m2   = state.fluid(2).m;
eta1 = state.fluid(1).eta;
eta2 = state.fluid(2).eta;
q1   = state.fluid(1).q;
q2   = state.fluid(2).q;
w1   = state.fluid(1).w;
w2   = state.fluid(2).w;
u1   = state.fluid(1).u;
u2   = state.fluid(2).u;
v1   = state.fluid(1).v;
v2   = state.fluid(2).v;
tke1 = state.fluid(1).tke;
tke2 = state.fluid(2).tke;
Tw11 = state.fluid(1).Tw(1);
Tw21 = state.fluid(2).Tw(1);
vareta1 = state.fluid(1).vareta;
vareta2 = state.fluid(2).vareta;
varq1   = state.fluid(1).varq;
varq2   = state.fluid(2).varq;
constants = settings.constants;
gravity = constants.phys.gravity;
bentraint = constants.param.bentraint;
bentrainq = constants.param.bentrainq;
bentrainw = constants.param.bentrainw;
bentrainu = constants.param.bentrainu;
bdetraint = constants.param.bdetraint;
bdetrainq = constants.param.bdetrainq;
bdetrainw = constants.param.bdetrainw;
bdetrainu = constants.param.bdetrainu;

% Remap mass to w levels
m1bar = weight_to_w(grid,m1);
m2bar = weight_to_w(grid,m2);

% Mean over the two fluids
rho    = m1 + m2;
rhobar = m1bar + m2bar;

% Fields needed to implement certain approximations
if switches.b | switches.e
   eta    = (m1bar.*eta1 + m2bar.*eta2)./rhobar;
   q      = (m1bar.*q1   + m2bar.*q2  )./rhobar;
   w      = (m1bar.*w1   + m2bar.*w2  )./rhobar;
   u      = (m1.*u1      + m2.*u2     )./rho;
   v      = (m1.*v1      + m2.*v2     )./rho;
end

% ------

% Surface pressure
psurf = grid.extrapb1*p(1) + grid.extrapb2*p(2);

% Vertical pressure gradient
dpdz(2:nz)   = (p(2:nz) - p(1:nz-1))./dzw(2:nz);
dpdz(1) = dpdz(2);
dpdz(nzp) = dpdz(nz);
dpdzbar = abovep.*dpdz(2:nzp) + belowp.*dpdz(1:nz);
dpdz(1) = 0;
dpdz(nzp) = 0;

% ------

% Diagnose densities, and compute residuals in equations of state
% and terms needed in linearization
if switches.a
    % Approximation a: Use mean fluid and fluid 2 thermodynamic equations
    eos = find_eos_a(grid, state, constants);
else
    % Use fluid 1 and fluid 2 thermodynamic equations
    %eos = find_eos(grid, state, constants);
    %work.nsq1 = eos.nsq1;
    %work.nsq2 = eos.nsq2;
    eos = find_eos_sg(grid, state, constants);
    % Copy for testing
%     eos.rho_deriv_eta1 = eosx.rho_deriv_eta1;
%     eos.rho_deriv_eta2 = eosx.rho_deriv_eta2;
%     eos.rho_deriv_q1 = eosx.rho_deriv_q1;
%     eos.rho_deriv_q2 = eosx.rho_deriv_q2;
%     eos.cldfrac1 = eosx.cldfrac1;
%     eos.cldfrac2 = eosx.cldfrac2;
%     eos.ql1 = eosx.ql1;
%     eos.ql2 = eosx.ql2;
end

% ------

% Factors needed for entrainment effects of vertical diffusion

% sigma2 * rho (normalized to make sure sum of sigmas = 1)
nearly_one = 1./(eos.sigma1 + eos.sigma2);
sig2rho = nearly_one.*eos.sigma2.*rho;

% And its weighted vertical average
sig2rhobar = weight_to_w(grid,sig2rho);
sig1rhobar = rhobar - sig2rhobar;

% d/dz of rho
ddz_rho(1) = 0;
ddz_rho(2:nz) = (rho(2:nz) - rho(1:nz-1))./dzw(2:nz);
ddz_rho(nzp) = 0;

% d/dz of sigma2*rho
ddz_sigrho(1) = 0;
ddz_sigrho(2:nz) = (sig2rho(2:nz) - sig2rho(1:nz-1))./dzw(2:nz);
ddz_sigrho(nzp) = 0;

% ------

% Determine surface forcing and the partitioning of surface fluxes
force = set_forcing(settings.forcing, t);
surface_flux = find_surface_flux(state,grid,eos,force,psurf,dpdzbar,constants);

% Calculate certain scales - these are non-local quantities, so difficult
% to linearize, but slowly evolving. (Should be OK to calculate forward in
% time (except first time step of the ARM case!!)
scales = find_scales(grid,settings,state,eos,surface_flux,constants,switches);

% ------
            
% Set diffusion coefficients

% Specified function of z and zstar ...
% [kdifft1,kdiffq1,kdiffw1,kdifft2,kdiffq2,kdiffw2,kdiffu1,kdiffu2,kdifftke1,kdifftke2] = ...
   % set_diffusion(grid.zp,grid.zw,scales.zstar,scales.wstar,scales.ustar_by_wstar,...
                 % constants.phys.k,constants.param.zrough,switches.d);
[kdifft1x,kdiffq1x,kdiffw1x,kdifft2x,kdiffq2x,kdiffw2x,kdiffu1x,kdiffu2x,kdifftke1x,kdifftke2x] = ...
    set_diffusion(grid.zp,grid.zw,scales.zstar,scales.wstar,scales.ustar_by_wstar,...
                  constants.phys.k,constants.param.zrough,switches.d);
                    
% ... or diagnosed from tke and turbulent length scales
[kdifft1,kdiffq1,kdiffw1,kdifft2,kdiffq2,kdiffw2,kdiffu1,kdiffu2,kdifftke1,kdifftke2] = ...
             set_diffusion_2(grid,scales.L_turb1,scales.L_turb2, ...
                             tke1,tke2,switches.d);
% [kdifft1,kdiffq1,kdiffw1,kdifft2,kdiffq2,kdiffw2,kdiffu1,kdiffu2,kdifftke1,kdifftke2] = ...
             % set_diffusion_2(grid,scales.L_turb1,scales.L_turb2, ...
                             % tke1.*m1./(m1+m2),tke2.*m2./(m1+m2),switches.d);
% [kdifft1,kdiffq1,kdiffw1,kdifft2,kdiffq2,kdiffw2,kdiffu1,kdiffu2,kdifftke1,kdifftke2] = ...
             % set_diffusion_2(grid,scales.L_turb1,scales.L_turb2, ...
                             % tke1.*m1./(m1+m2)+tke2.*m2./(m1+m2),tke1.*m1./(m1+m2)+tke2.*m2./(m1+m2),switches.d);

% For testing...
% If old kdiff has been supplied then use it              
%if old_diff.flag
%    disp('** Using old Kdiffs **')
%    kdifft1 = old_diff.kdifft1;
%    kdifft2 = old_diff.kdifft2;
%    kdiffq1 = old_diff.kdiffq1;
%    kdiffq2 = old_diff.kdiffq2;
%    kdiffw1 = old_diff.kdiffw1;
%    kdiffw2 = old_diff.kdiffw2;
%    kdiffu1 = old_diff.kdiffu1;
%    kdiffu2 = old_diff.kdiffu2;
%    kdifftke1 = old_diff.kdifftke1;
%    kdifftke2 = old_diff.kdifftke2;
%end

% ------

% Advective mass flux
[F1, dF1dma, dF1dmb, dF1dw] = mass_flux(grid,m1,w1);
[F2, dF2dma, dF2dmb, dF2dw] = mass_flux(grid,m2,w2);

% Include surface moisture fluxes
F1(1) = F1(1) + surface_flux.q1;
F2(1) = F2(1) + surface_flux.q2;

% Mass tendencies
tend.fluid(1).m.transport = - (F1(2:nzp) - F1(1:nz))./dzp;
tend.fluid(2).m.transport = - (F2(2:nzp) - F2(1:nz))./dzp;

% ------

% Interpolate mass fluxes to p levels
F1bar = abovep.*F1(2:nzp) + belowp.*F1(1:nz);
F2bar = abovep.*F2(2:nzp) + belowp.*F2(1:nz);

% ------

% Advective entropy flux
if switches.e
    [Feta1,dFeta1detaa,dFeta1detab,deta1udz,eta1ubar] = tracer_flux(grid,F1bar,eta );
else
    [Feta1,dFeta1detaa,dFeta1detab,deta1udz,eta1ubar] = tracer_flux(grid,F1bar,eta1);
end
[Feta2,dFeta2detaa,dFeta2detab,deta2udz,eta2ubar] = tracer_flux(grid,F2bar,eta2);

% Transport tendencies of mass times entropy
% Boundary advective fluxes are assumed zero
tend.fluid(1).meta.transport(1:nz ) = - Feta1./dzw(1:nz);
tend.fluid(1).meta.transport(nzp) = 0;
tend.fluid(1).meta.transport(2:nzp) = tend.fluid(1).meta.transport(2:nzp) ...
                                      + Feta1./dzw(2:nzp);
tend.fluid(2).meta.transport(1:nz ) = - Feta2./dzw(1:nz);
tend.fluid(2).meta.transport(nzp) = 0;
tend.fluid(2).meta.transport(2:nzp) = tend.fluid(2).meta.transport(2:nzp) ...
                                      + Feta2./dzw(2:nzp);

% Diffusive entropy flux
[Deta2, dDeta2detaa, dDeta2detab, dDeta2dm ] = diff_flux( grid, kdifft2, eta2 , m2, surface_flux.eta2);
if switches.b
    [Deta1, dDeta1detaa, dDeta1detab, dDeta1dm ] = diff_flux( grid, kdifft1, eta  , rho, surface_flux.eta);
    Deta1 = Deta1 - Deta2;
    dDeta1detaa = dDeta1detaa.*m1bar(2:nzp)./rhobar(2:nzp);
    dDeta1detab = dDeta1detab.*m1bar(1:nz )./rhobar(1:nz );
    % dDeta1dm    = dDeta1dm;
else
    [Deta1, dDeta1detaa, dDeta1detab, dDeta1dm ] = diff_flux( grid, kdifft1, eta1 , m1, surface_flux.eta1);
end
% Include contribution from surface flux to first interior level
Deta1(1) = Deta1(1) + surface_flux.L1eta1;
Deta2(1) = Deta2(1) + surface_flux.L1eta2;
% disp(['Deta1 orig = ',num2str(Deta1(1:4))])

% Include buoyancy correlation term
rr = 0.0;
covaretaq1 = rr*sqrt(vareta1.*varq1);
covaretaq2 = rr*sqrt(vareta2.*varq2);
t_scale1 = 6*scales.L_turb1./sqrt(tke1);
t_scale2 = 6*scales.L_turb2./sqrt(tke2);
Bterm = eos.rho_deriv_eta1.*vareta1 + eos.rho_deriv_q1.*covaretaq1;
Deta1bc = - gravity*t_scale1.*eos.sigma1.*(abovep.*Bterm(2:nzp) + belowp.*Bterm(1:nz));
Bterm = eos.rho_deriv_eta2.*vareta2 + eos.rho_deriv_q2.*covaretaq2;
Deta2bc = - gravity*t_scale2.*eos.sigma2.*(abovep.*Bterm(2:nzp) + belowp.*Bterm(1:nz));
% disp(['Deta1 new  = ',num2str(Deta1(1:4))])
% disp(['L_turb1        = ',num2str(scales.L_turb1(1:4))])
% disp(['U1             = ',num2str(sqrt(tke1(1:4)))])
% disp(['t_scale1       = ',num2str(t_scale1(1:4))])
% disp(['K              = ',num2str(kdifft1(1:4))])
% disp(['rho_deriv_eta1 = ',num2str(eos.rho_deriv_eta1(1:4))])
% disp(['vareta1        = ',num2str(vareta1(1:4))])


% Diffusive tendencies of mass times entropy
tend.fluid(1).meta.diffuse(1:nz ) = - Deta1./dzw(1:nz);
tend.fluid(1).meta.diffuse(nzp) = 0;
tend.fluid(1).meta.diffuse(2:nzp) = tend.fluid(1).meta.diffuse(2:nzp) ...
                                    + Deta1./dzw(2:nzp);
tend.fluid(2).meta.diffuse(1:nz ) = - Deta2./dzw(1:nz);
tend.fluid(2).meta.diffuse(nzp) = 0;
tend.fluid(2).meta.diffuse(2:nzp) = tend.fluid(2).meta.diffuse(2:nzp) ...
                                    + Deta2./dzw(2:nzp);

% Include surface entropy fluxes,
tend.fluid(1).meta.diffuse(1) = tend.fluid(1).meta.diffuse(1) + surface_flux.eta1/dzw(1);
tend.fluid(2).meta.diffuse(1) = tend.fluid(2).meta.diffuse(1) + surface_flux.eta2/dzw(1);

% Entrainment/detrainment associated with vertical diffusion and gradients
% of sigma
if switches.b
    disp('*** Need to formulate diffent term for switch b ***')
    pause
end
% Mean K grad eta over the two fluids
meanG = -(Deta1 + Deta2)./rho;
% Vertical average
meanGbar(1:nz) = grid.abovew(1:nz).*meanG;
meanGbar(nzp) = 0;
meanGbar(2:nzp) = meanGbar(2:nzp) + grid.beloww(2:nzp).*meanG;
% Total correction
corrde = meanGbar.*(ddz_sigrho - (m2bar./rhobar).*ddz_rho);
tend.fluid(1).meta.diffent =   corrde;
tend.fluid(2).meta.diffent = - corrde;

% ------

% Advective water flux
if switches.e
    [Fq1,dFq1dqa,dFq1dqb,dq1udz,q1ubar] = tracer_flux(grid,F1bar,q );
else
    [Fq1,dFq1dqa,dFq1dqb,dq1udz,q1ubar] = tracer_flux(grid,F1bar,q1);
end
[Fq2,dFq2dqa,dFq2dqb,dq2udz,q2ubar] = tracer_flux(grid,F2bar,q2);

% Transport Tendencies of mass times q
tend.fluid(1).mq.transport(1:nz ) = - Fq1./dzw(1:nz);
tend.fluid(1).mq.transport(nzp) = 0;
tend.fluid(1).mq.transport(2:nzp) = tend.fluid(1).mq.transport(2:nzp) ...
                        + Fq1./dzw(2:nzp);
tend.fluid(2).mq.transport(1:nz ) = - Fq2./dzw(1:nz);
tend.fluid(2).mq.transport(nzp) = 0;
tend.fluid(2).mq.transport(2:nzp) = tend.fluid(2).mq.transport(2:nzp) ...
                        + Fq2./dzw(2:nzp);
% Include surface moisture fluxes,
tend.fluid(1).mq.transport(1) = tend.fluid(1).mq.transport(1) + surface_flux.q1/dzw(1);
tend.fluid(2).mq.transport(1) = tend.fluid(2).mq.transport(1) + surface_flux.q2/dzw(1);

% Diffusive water flux
[Dq2, dDq2dqa, dDq2dqb, dDq2dm ] = diff_flux( grid, kdiffq2, q2 , m2, surface_flux.q2);
if switches.b
    [Dq1, dDq1dqa, dDq1dqb, dDq1dm ] = diff_flux( grid, kdiffq1, q  , rho, surface_flux.q);
    Dq1 = Dq1 - Dq2;
    dDq1dqa = dDq1dqa.*m1bar(2:nzp)./rhobar(2:nzp);
    dDq1dqb = dDq1dqb.*m1bar(1:nz )./rhobar(1:nz );
    % dDq1dm  = dDq1dm;
else
    [Dq1, dDq1dqa, dDq1dqb, dDq1dm ] = diff_flux( grid, kdiffq1, q1 , m1, surface_flux.q1);
end
% Include contribution from surface flux to first interior level
Dq1(1) = Dq1(1) + surface_flux.L1q1;
Dq2(1) = Dq2(1) + surface_flux.L1q2;

% Include buoyancy correlation term
% Correlations and timescales are computed above
% rr = 0.0;
% covaretaq1 = rr*sqrt(vareta1.*varq1);
% covaretaq2 = rr*sqrt(vareta2.*varq2);
% t_scale1 = 6*scales.L_turb1./sqrt(tke1);
% t_scale2 = 6*scales.L_turb2./sqrt(tke2);
% Bterm = eos.rho_deriv_eta1.*covaretaq1 + eos.rho_deriv_q1.*varq1;
% Dq1 = Dq1 - gravity*t_scale1.*eos.sigma1.*(abovep.*Bterm(2:nzp) + belowp.*Bterm(1:nz));
% Bterm = eos.rho_deriv_eta2.*covaretaq2 + eos.rho_deriv_q2.*varq2;
% Dq2 = Dq2 - gravity*t_scale2.*eos.sigma2.*(abovep.*Bterm(2:nzp) + belowp.*Bterm(1:nz));

% Tendencies of mass times q
tend.fluid(1).mq.diffuse(1:nz ) = - Dq1./dzw(1:nz);
tend.fluid(1).mq.diffuse(nzp) = 0;
tend.fluid(1).mq.diffuse(2:nzp) = tend.fluid(1).mq.diffuse(2:nzp) ...
                        + Dq1./dzw(2:nzp);
tend.fluid(2).mq.diffuse(1:nz ) = - Dq2./dzw(1:nz);
tend.fluid(2).mq.diffuse(nzp) = 0;
tend.fluid(2).mq.diffuse(2:nzp) = tend.fluid(2).mq.diffuse(2:nzp) ...
                        + Dq2./dzw(2:nzp);

% Entrainment/detrainment associated with vertical diffusion and gradients
% of sigma
if switches.b
    disp('*** Need to formulate diffent term for switch b ***')
    pause
end
% Mean K grad eta over the two fluids
meanG = -(Dq1 + Dq2)./rho;
% Vertical average
meanGbar(1:nz) = grid.abovew(1:nz).*meanG;
meanGbar(nzp) = 0;
meanGbar(2:nzp) = meanGbar(2:nzp) + grid.beloww(2:nzp).*meanG;
% Total correction
corrde = meanGbar.*(ddz_sigrho - (m2bar./rhobar).*ddz_rho);
tend.fluid(1).mq.diffent =   corrde;
tend.fluid(2).mq.diffent = - corrde;

% ------

% Non-hydrostatic pressure gradient terms
nhpg1 = dpdz./eos.rhow1 + gravity;
nhpg2 = dpdz./eos.rhow2 + gravity;

% Diffusive w flux
[Dw2, dDw2dwa, dDw2dwb, dDw2dm ] = diff_flux( grid, kdiffw2, w2 , m2, 0);
if switches.b
    [Dw1, dDw1dwa, dDw1dwb, dDw1dm ] = diff_flux( grid, kdiffw1, w  , rho, 0);
    Dw1 = Dw1 - Dw2;
    dDw1dwa = dDw1dwa.*m1bar(2:nzp)./rhobar(2:nzp);
    dDw1dwb = dDw1dwb.*m1bar(1:nz )./rhobar(1:nz );
else
    [Dw1, dDw1dwa, dDw1dwb, dDw1dm ] = diff_flux( grid, kdiffw1, w1 , m1, 0);
end

% Entrainment/detrainment associated with vertical diffusion and gradients
% of sigma
if switches.b
    disp('*** Need to formulate diffent term for switch b ***')
    pause
end
% Mean K grad eta over the two fluids
meanG = -(Dw1 + Dw2)./rho;
% Vertical average
meanGbar(1:nz) = grid.abovew(1:nz).*meanG;
meanGbar(nzp) = 0;
meanGbar(2:nzp) = meanGbar(2:nzp) + grid.beloww(2:nzp).*meanG;
% Total correction
corrde = meanGbar.*(ddz_sigrho - (m2bar./rhobar).*ddz_rho);

% Advective w flux
if switches.e
    [Fw1,dFw1dwa,dFw1dwb,dw1udz,w1ubar] = tracer_flux(grid,F1bar,w );
else
    [Fw1,dFw1dwa,dFw1dwb,dw1udz,w1ubar] = tracer_flux(grid,F1bar,w1);
end
[Fw2,dFw2dwa,dFw2dwb,dw2udz,w2ubar] = tracer_flux(grid,F2bar,w2);

% Transport Tendencies of mass times w
tend.fluid(1).mw.transport(1:nz ) = - Fw1./dzw(1:nz);
tend.fluid(1).mw.transport(nzp) = 0;
tend.fluid(1).mw.transport(2:nzp) = tend.fluid(1).mw.transport(2:nzp) ...
                        + Fw1./dzw(2:nzp);
tend.fluid(2).mw.transport(1:nz ) = - Fw2./dzw(1:nz);
tend.fluid(2).mw.transport(nzp) = 0;
tend.fluid(2).mw.transport(2:nzp) = tend.fluid(2).mw.transport(2:nzp) ...
                        + Fw2./dzw(2:nzp);

% Length scale for pressure drag
zdrag = scales.zstar;
% Pressure drag
[ drag, ddragdw1, ddragdw2 ] = find_drag( m2bar, w1, w2, zdrag );

% Other tendencies of mass times w
tend.fluid(1).mw.drag = drag;
tend.fluid(1).mw.diffuse(1) = 0;
tend.fluid(1).mw.diffuse(2:nz) = (Dw1(1:nz-1) - Dw1(2:nz))./dzw(2:nz);
tend.fluid(1).mw.diffuse(nzp) = 0;
tend.fluid(1).mw.pgterm = -m1bar.*nhpg1;
tend.fluid(1).mw.pgterm(1) = 0;
tend.fluid(1).mw.pgterm(nzp) = 0;
tend.fluid(2).mw.drag = - drag;
tend.fluid(2).mw.diffuse(1) = 0;
tend.fluid(2).mw.diffuse(2:nz) = (Dw2(1:nz-1) - Dw2(2:nz))./dzw(2:nz);
tend.fluid(2).mw.diffuse(nzp) = 0;
tend.fluid(1).mw.diffent =   corrde;
tend.fluid(2).mw.diffent = - corrde;
tend.fluid(2).mw.pgterm = -m2bar.*nhpg2;
tend.fluid(2).mw.pgterm(1) = 0;
tend.fluid(2).mw.pgterm(nzp) = 0;

% ------

% Coriolis terms
tend.fluid(1).mu.coriolis =   constants.phys.coriolis*m1.*(state.fluid(1).v - force.vg);
tend.fluid(2).mu.coriolis =   constants.phys.coriolis*m2.*(state.fluid(2).v - force.vg);
tend.fluid(1).mv.coriolis = - constants.phys.coriolis*m1.*(state.fluid(1).u - force.ug);
tend.fluid(2).mv.coriolis = - constants.phys.coriolis*m2.*(state.fluid(2).u - force.ug);

% Advective u and v fluxes
if switches.e
    [Fu1,dFu1dua,dFu1dub,du1udz,u1ubar] = tracer_flux_p(grid,F1,u );
    [Fv1,dFv1dva,dFv1dvb,dv1udz,v1ubar] = tracer_flux_p(grid,F1,v );
else
    [Fu1,dFu1dua,dFu1dub,du1udz,u1ubar] = tracer_flux_p(grid,F1,u1);
    [Fv1,dFv1dva,dFv1dvb,dv1udz,v1ubar] = tracer_flux_p(grid,F1,v1);
end
[Fu2,dFu2dua,dFu2dub,du2udz,u2ubar] = tracer_flux_p(grid,F2,u2);
[Fv2,dFv2dva,dFv2dvb,dv2udz,v2ubar] = tracer_flux_p(grid,F2,v2);
tend.fluid(1).mu.transport = -(Fu1(2:nzp) - Fu1(1:nz))./dzp;
tend.fluid(1).mv.transport = -(Fv1(2:nzp) - Fv1(1:nz))./dzp;
tend.fluid(2).mu.transport = -(Fu2(2:nzp) - Fu2(1:nz))./dzp;
tend.fluid(2).mv.transport = -(Fv2(2:nzp) - Fv2(1:nz))./dzp;

% Diffusive u and v fluxes
speed = sqrt(u1(1)*u1(1) + v1(1)*v1(1));
[Du2, dDu2dua, dDu2dub ] = diff_flux_u( grid, kdiffu2, u2 , m2bar, speed, constants.phys.k, constants.param.zrough);
[Dv2, dDv2dva, dDv2dvb ] = diff_flux_u( grid, kdiffu2, v2 , m2bar, speed, constants.phys.k, constants.param.zrough);
if switches.b
    [Du1, dDu1dua, dDu1dub ] = diff_flux_u( grid, kdiffu1, u  , rhobar, speed, constants.phys.k, constants.param.zrough);
    [Dv1, dDv1dva, dDv1dvb ] = diff_flux_u( grid, kdiffu1, v  , rhobar, speed, constants.phys.k, constants.param.zrough);
    Du1 = Du1 - Du2;
    dDu1dua = dDu1dua - dDu2dua;
    dDu1dub = dDu1dub - dDu2dub;
    Dv1 = Dv1 - Dv2;
    dDv1dva = dDv1dva - dDv2dva;
    dDv1dvb = dDv1dvb - dDv2dvb;
else
    [Du1, dDu1dua, dDu1dub ] = diff_flux_u( grid, kdiffu1, u1 , m1bar, speed, constants.phys.k, constants.param.zrough);
    [Dv1, dDv1dva, dDv1dvb ] = diff_flux_u( grid, kdiffu1, v1 , m1bar, speed, constants.phys.k, constants.param.zrough);
end
tend.fluid(1).mu.diffuse = -(Du1(2:nzp) - Du1(1:nz))./dzp;
tend.fluid(1).mv.diffuse = -(Dv1(2:nzp) - Dv1(1:nz))./dzp;
tend.fluid(2).mu.diffuse = -(Du2(2:nzp) - Du2(1:nz))./dzp;
tend.fluid(2).mv.diffuse = -(Dv2(2:nzp) - Dv2(1:nz))./dzp;

% Entrainment/detrainment associated with vertical diffusion and gradients
% of sigma
if switches.b
    disp('*** Need to formulate diffent term for switch b ***')
    pause
end
% Mean K grad eta over the two fluids
meanG = -(Du1 + Du2)./rhobar;
s1rG = sig1rhobar.*meanG;
s2rG = sig2rhobar.*meanG;
corrde = nearly_one ...
       .*( eos.sigma1.*(s2rG(2:nzp) - s2rG(1:nz)) ...
         - eos.sigma2.*(s1rG(2:nzp) - s1rG(1:nz)) )./dzp;
tend.fluid(1).mu.diffent =   corrde;
tend.fluid(2).mu.diffent = - corrde;

% Entrainment/detrainment associated with vertical diffusion and gradients
% of sigma
if switches.b
    disp('*** Need to formulate diffent term for switch b ***')
    pause
end
% Mean K grad eta over the two fluids
meanG = -(Dv1 + Dv2)./rhobar;
s1rG = sig1rhobar.*meanG;
s2rG = sig2rhobar.*meanG;
corrde = nearly_one ...
       .*( eos.sigma1.*(s2rG(2:nzp) - s2rG(1:nz)) ...
         - eos.sigma2.*(s1rG(2:nzp) - s1rG(1:nz)) )./dzp;
tend.fluid(1).mv.diffent =   corrde;
tend.fluid(2).mv.diffent = - corrde;

% Neglect horizontal inter-fluid PG terms for now
% Also neglect interfluid forces due to KH instability (??)

% ------

% Advective tke fluxes
[Ftke1,dFtke1dtkea,dFtke1dtkeb,dtke1udz,tke1ubar] = tracer_flux_p(grid,F1,tke1);
[Ftke2,dFtke2dtkea,dFtke2dtkeb,dtke2udz,tke2ubar] = tracer_flux_p(grid,F2,tke2);
tend.fluid(1).mtke.transport = -(Ftke1(2:nzp) - Ftke1(1:nz))./dzp;
tend.fluid(2).mtke.transport = -(Ftke2(2:nzp) - Ftke2(1:nz))./dzp;

% Diffusive tke fluxes
[Dtke1, dDtke1dtkea, dDtke1dtkeb ] = diff_flux_tke( grid, kdifftke1, tke1 , m1bar);
[Dtke2, dDtke2dtkea, dDtke2dtkeb ] = diff_flux_tke( grid, kdifftke2, tke2 , m2bar);
tend.fluid(1).mtke.diffuse = -(Dtke1(2:nzp) - Dtke1(1:nz))./dzp;
tend.fluid(2).mtke.diffuse = -(Dtke2(2:nzp) - Dtke2(1:nz))./dzp;

% Entrainment/detrainment associated with vertical diffusion and gradients
% of sigma
if switches.b
    disp('*** Need to formulate diffent term for switch b ***')
    pause
end
% Mean K grad eta over the two fluids
meanG = -(Dtke1 + Dtke2)./rhobar;
s1rG = sig1rhobar.*meanG;
s2rG = sig2rhobar.*meanG;
corrde = nearly_one ...
       .*( eos.sigma1.*(s2rG(2:nzp) - s2rG(1:nz)) ...
         - eos.sigma2.*(s1rG(2:nzp) - s1rG(1:nz)) )./dzp;
tend.fluid(1).mtke.diffent =   corrde;
tend.fluid(2).mtke.diffent = - corrde;

% Shear generation
% Consistent with free slip at top boundary
corrdew       = (w2 - w1).*tend.fluid(1).mw.diffent;
corrdew       = grid.aboves.*corrdew(2:nzp) + grid.belows.*corrdew(1:nz);
corrdeu       = (u2 - u1).*tend.fluid(1).mu.diffent;
corrdev       = (v2 - v1).*tend.fluid(1).mv.diffent;
dudz(1)       =  u1(1)/dzw(1);
dudz(2:nz)    = (u1(2:nz) - u1(1:nz-1))./dzw(2:nz);
dudz(nzp)     =  0;
dvdz(1)       =  v1(1)/dzw(1);
dvdz(2:nz)    = (v1(2:nz) - v1(1:nz-1))./dzw(2:nz);
dvdz(nzp)     =  0;
kesinkw       = - Dw1.*(w1(2:nzp) - w1(1:nz))./dzp ...
                + eos.sigma1.*corrdew;
Kdudzsq       = kdiffu1.*(dudz.*dudz + dvdz.*dvdz);
Kdudzsq(1)    = - (u1(1)*Du1(1) + v1(1)*Dv1(1))/(m1bar(1)*dzw(1));
kesinkuv      = m1.*(grid.aboves.*Kdudzsq(2:nzp) + grid.belows.*Kdudzsq(1:nz)) ...
              + eos.sigma1.*(corrdeu + corrdev);
tend.fluid(1).mtke.shear = kesinkw + kesinkuv;
dudz(1)       =  u2(1)/dzw(1);
dudz(2:nz)    = (u2(2:nz) - u2(1:nz-1))./dzw(2:nz);
dudz(nzp)     =  0;
dvdz(1)       =  v2(1)/dzw(1);
dvdz(2:nz)    = (v2(2:nz) - v2(1:nz-1))./dzw(2:nz);
dvdz(nzp)     =  0;
kesinkw       = - Dw2.*(w2(2:nzp) - w2(1:nz))./dzp ...
                + eos.sigma2.*corrdew;
Kdudzsq       = kdiffu2.*(dudz.*dudz + dvdz.*dvdz);
Kdudzsq(1)    = - (u2(1)*Du2(1) + v2(1)*Dv2(1))/(m2bar(1)*dzw(1));
kesinkuv      = m2.*(grid.aboves.*Kdudzsq(2:nzp) + grid.belows.*Kdudzsq(1:nz)) ...
              + eos.sigma2.*(corrdeu + corrdev);
tend.fluid(2).mtke.shear = kesinkw + kesinkuv;

% Buoyancy generation
% Limit to prevent tke sink when there isn't enough tke
bflux1 = dpdzbar.*(eos.drdetap1.*Deta1 + eos.drdqp1.*Dq1);
tend.fluid(1).mtke.bflux = max(bflux1,m1.*(constants.param.tke_min - tke1)/dt);
bflux2 = dpdzbar.*(eos.drdetap2.*Deta2 + eos.drdqp2.*Dq2);
tend.fluid(2).mtke.bflux = max(bflux2,m2.*(constants.param.tke_min - tke2)/dt);

% Drag term
% It is not obvious how to partition this between fluid 1 and fluid 2
% so just divide it evenly
kesinkw = - (w1 - w2).*drag;
tkesrc = grid.aboves.*kesinkw(2:nzp) + grid.belows.*kesinkw(1:nz);
tend.fluid(1).mtke.drag = 0.5*tkesrc;
tend.fluid(2).mtke.drag = 0.5*tkesrc;

% Dissipation term
tend.fluid(1).mtke.dissn = - m1.*tke1./scales.T_turb1;
tend.fluid(2).mtke.dissn = - m2.*tke2./scales.T_turb2;

% ------

% Include effect of buoyancy flux on internal energy
esource = weight_to_w(grid,tend.fluid(1).mtke.bflux./m1);
tend.fluid(1).meta.bflux = -m1bar.*esource./state.fluid(1).Tw;
esource = weight_to_w(grid,tend.fluid(2).mtke.bflux./m2);
tend.fluid(2).meta.bflux = -m2bar.*esource./state.fluid(2).Tw;

% Return dissipated TKE as entropy
esource = weight_to_w(grid,tend.fluid(1).mtke.dissn./m1);
tend.fluid(1).meta.dissn = -m1bar.*esource./state.fluid(1).Tw;
esource = weight_to_w(grid,tend.fluid(2).mtke.dissn./m2);
tend.fluid(2).meta.dissn = -m2bar.*esource./state.fluid(2).Tw;

% Note entropy source due to diffusion is not yet included, but
% it is estimated to be small

% ------

% Turbulence time scales at w-levels
T_turb1_bar = weight_to_w(grid,scales.T_turb1);
T_turb2_bar = weight_to_w(grid,scales.T_turb2);

%disp(['T turb 1 = ',num2str(T_turb1_bar(1:4))])

% ------

% eta variance
% Neglect transport terms (MYNN level 2.5)
% ?? Would it be better to store variances at p levels ??

% Downgradient flux source terms ...
deta1dz = (eta1(2:nzp) - eta1(1:nz))./grid.dzp;
deta2dz = (eta2(2:nzp) - eta2(1:nz))./grid.dzp;
dg1 = -Deta1.*deta1dz;
dg2 = -Deta2.*deta2dz;
% ... mapped to w levels
tend.fluid(1).mvareta.diffuse = weight_to_w(grid,dg1);
tend.fluid(2).mvareta.diffuse = weight_to_w(grid,dg2);

% Buoyancy correlation terms
% disp('** No Buoyancy correl in eta variance budget **')
% Bterm = eos.rho_deriv_eta1.*vareta1 + eos.rho_deriv_q1.*covaretaq1;
% tend.fluid(1).mvareta.diffuse = tend.fluid(1).mvareta.diffuse ...
%     + weight_to_w(grid,gravity*t_scale1.*eos.sigma1.*deta1dz).*Bterm;
% Bterm = eos.rho_deriv_eta2.*vareta2 + eos.rho_deriv_q2.*covaretaq2;
% tend.fluid(2).mvareta.diffuse = tend.fluid(2).mvareta.diffuse ...
%     + weight_to_w(grid,gravity*t_scale2.*eos.sigma2.*deta2dz).*Bterm;

% `diffent' source terms
corrde = (eta2 - eta1).*tend.fluid(1).meta.diffent;
tend.fluid(1).mvareta.diffent = m1bar.*corrde./rhobar;
tend.fluid(2).mvareta.diffent = m2bar.*corrde./rhobar;

% Dissipation terms
tend.fluid(1).mvareta.dissn = - m1bar.*vareta1./T_turb1_bar;
tend.fluid(2).mvareta.dissn = - m2bar.*vareta2./T_turb2_bar;

% ------

% q variance
% Neglect transport terms (MYNN level 2.5)

% Downgradient flux source terms ...
dq1dz = (q1(2:nzp) - q1(1:nz))./grid.dzp;
dq2dz = (q2(2:nzp) - q2(1:nz))./grid.dzp;
dg1 = -Dq1.*dq1dz;
dg2 = -Dq2.*dq2dz;
% ... mapped to w levels
tend.fluid(1).mvarq.diffuse = weight_to_w(grid,dg1);
tend.fluid(2).mvarq.diffuse = weight_to_w(grid,dg2);

% `diffent' source terms
corrde = (q2 - q1).*tend.fluid(1).mq.diffent;
tend.fluid(1).mvarq.diffent = m1bar.*corrde./rhobar;
tend.fluid(2).mvarq.diffent = m2bar.*corrde./rhobar;

% Dissipation terms
tend.fluid(1).mvarq.dissn = - m1bar.*varq1./T_turb1_bar;
tend.fluid(2).mvarq.dissn = - m2bar.*varq2./T_turb2_bar;

% ------

% Now include entrainment/detrainment contributions

% Diagnose buoyancy (needed for diagnostics and maybe for relabelling)
sigma1 = m1./eos.rho1;
sigma1bar(2:nz) = grid.abovew(2:nz).*sigma1(2:nz) + grid.beloww(2:nz).*sigma1(1:nz-1);
sigma1bar(1)   = sigma1bar(2);
sigma1bar(nzp) = sigma1bar(nz);

% rho_mean = weight_to_w(grid, m1 + m2);
% buoy = gravity*sigma1bar.*(eos.rhow1 - rho_mean)./rho_mean;

buoy = gravity*sigma1bar.*(eos.rhow1 - eos.rhow2)./eos.rhow2;
%buoyx = gravity*sigma1bar.*(eosx.rhow1 - eosx.rhow2)./eosx.rhow2;
% Calculate resolved buoyancy flux while we're here
buoy_flux_res = F2.*buoy./sigma1bar;

% Set entrainment/detrainment coefficients
% Trial version
relabel  = set_entrain_trial     (grid, state, buoy, eos, scales, kdiffw2, constants, dt);
%relabelx = set_entrain_trial_save(grid, state, buoy, eos, scales, kdiffw2, constants, dt);
%compare_relabel


% Unpack entrainment variables
M12 = relabel.M12;
M21 = relabel.M21;
M12bar = relabel.M12bar;
M21bar = relabel.M21bar;

% Mass tendencies
tend.fluid(1).m.relabel = M12 - M21;
tend.fluid(2).m.relabel = M21 - M12;

% Mass times eta tendencies
tend.fluid(1).meta.relabel = M12bar.*relabel.etahat12 - M21bar.*relabel.etahat21;
tend.fluid(2).meta.relabel = M21bar.*relabel.etahat21 - M12bar.*relabel.etahat12;

% Mass times q tendencies
tend.fluid(1).mq.relabel = M12bar.*relabel.qhat12 - M21bar.*relabel.qhat21;
tend.fluid(2).mq.relabel = M21bar.*relabel.qhat21 - M12bar.*relabel.qhat12;

% Mass times w tendencies
tend.fluid(1).mw.relabel = (M12bar.*relabel.what12 - M21bar.*relabel.what21);
tend.fluid(2).mw.relabel = (M21bar.*relabel.what21 - M12bar.*relabel.what12);

% Mass times u and v tendencies
tend.fluid(1).mu.relabel = (M12.*relabel.uhat12 - M21.*relabel.uhat21);
tend.fluid(2).mu.relabel = (M21.*relabel.uhat21 - M12.*relabel.uhat12);
tend.fluid(1).mv.relabel = (M12.*relabel.vhat12 - M21.*relabel.vhat21);
tend.fluid(2).mv.relabel = (M21.*relabel.vhat21 - M12.*relabel.vhat12);

% Entrained and detrained values of TKE
% transferred tke
tkehat12 = tke2;
tkehat21 = tke1;
% % Resolved KE per unit mass
% wsq = w1.*w1;
% RKE1 = 0.5*(u1.*u1 + v1.*v1 + grid.aboves.*wsq(2:nzp) + grid.belows.*wsq(1:nz));
% wsq = w2.*w2;
% RKE2 = 0.5*(u2.*u2 + v2.*v2 + grid.aboves.*wsq(2:nzp) + grid.belows.*wsq(1:nz));
% % u dot transferred u
% wwhat = w1.*relabel.what12;
% u1uhat12 = u1.*relabel.uhat12 + v1.*relabel.vhat12 + grid.aboves.*wwhat(2:nzp) + grid.belows.*wwhat(1:nz);
% wwhat = w2.*relabel.what12;
% u2uhat12 = u2.*relabel.uhat12 + v2.*relabel.vhat12 + grid.aboves.*wwhat(2:nzp) + grid.belows.*wwhat(1:nz);
% wwhat = w2.*relabel.what21;
% u2uhat21 = u2.*relabel.uhat21 + v2.*relabel.vhat21 + grid.aboves.*wwhat(2:nzp) + grid.belows.*wwhat(1:nz);
% wwhat = w1.*relabel.what21;
% u1uhat21 = u1.*relabel.uhat21 + v1.*relabel.vhat21 + grid.aboves.*wwhat(2:nzp) + grid.belows.*wwhat(1:nz);
% % TKE tendencies due to entrainment/detrainment
% tend.fluid(1).mtke.relabel = M12.*(tkehat12 +   RKE1 + RKE2 - u1uhat12) ...
%                            - M21.*(tkehat21 + 2*RKE1        - u1uhat21);
% tend.fluid(2).mtke.relabel = M21.*(tkehat21 +   RKE2 + RKE1 - u2uhat21) ...
%                            - M12.*(tkehat12 + 2*RKE2        - u2uhat12);
dwsq = (relabel.what12 - w1).^2;
du12_1_sq = (relabel.uhat12 - u1).^2 ...
          + (relabel.vhat12 - v1).^2 ...
          + grid.aboves.*dwsq(2:nzp) + grid.belows.*dwsq(1:nz);
dwsq = (relabel.what21 - w1).^2;
du21_1_sq = (relabel.uhat21 - u1).^2 ...
          + (relabel.vhat21 - v1).^2 ...
          + grid.aboves.*dwsq(2:nzp) + grid.belows.*dwsq(1:nz);
dwsq = (relabel.what12 - w2).^2;
du12_2_sq = (relabel.uhat12 - u2).^2 ...
          + (relabel.vhat12 - v2).^2 ...
          + grid.aboves.*dwsq(2:nzp) + grid.belows.*dwsq(1:nz);
dwsq = (relabel.what21 - w2).^2;
du21_2_sq = (relabel.uhat21 - u2).^2 ...
          + (relabel.vhat21 - v2).^2 ...
          + grid.aboves.*dwsq(2:nzp) + grid.belows.*dwsq(1:nz);
% disp('*** Alternative TKE relabelling ***')
% TKE tendencies due to entrainment/detrainment
tend.fluid(1).mtke.relabel = M12.*(tkehat12 + 0.5*du12_1_sq) ...
                           - M21.*(tkehat21 + 0.5*du21_1_sq);
tend.fluid(2).mtke.relabel = M21.*(tkehat21 + 0.5*du21_2_sq) ...
                           - M12.*(tkehat12 + 0.5*du12_2_sq);

% Entrained and detrained values of eta variance
% Note these are quasi-advective form tendencies
% Resolved variance
Rvareta1 = eta1.*eta1;
Rvareta2 = eta2.*eta2;
% `Upwind' approximation for Transferred variances
% etavarhat12 = vareta1;
% etavarhat21 = vareta2;
% eta times transferred eta
%eta1etahat12 = eta1.*relabel.etahat12;
%eta2etahat12 = eta2.*relabel.etahat12;
%eta2etahat21 = eta2.*relabel.etahat21;
%eta1etahat21 = eta1.*relabel.etahat21;
% eta variance tendencies due to entrainment/detrainment
% tend.fluid(1).mvareta.relabel = M12bar.*(vareta2 - vareta1 ...
%                                        + Rvareta1 + Rvareta2 - 2*eta1etahat12) ...
%                               - M21bar.*(2*Rvareta1 - 2*eta1etahat21);
% tend.fluid(2).mvareta.relabel = M21bar.*(vareta1 - vareta2 ...
%                                        + Rvareta2 + Rvareta1 - 2*eta2etahat21) ...
%                               - M12bar.*(2*Rvareta2 - 2*eta2etahat12);
tend.fluid(1).mvareta.relabel = M12bar.*(vareta2 - vareta1 ...
                                      + (relabel.etahat12 - eta1).^2) ...
                              - M21bar.*(relabel.etahat21 - eta1).^2;
tend.fluid(2).mvareta.relabel = M21bar.*(vareta1 - vareta2 ...
                                      + (relabel.etahat21 - eta2).^2) ...
                              - M12bar.*(relabel.etahat12 - eta2).^2;
% disp(['M221,12   = ',num2str(M21bar(28)),' ',...
%                      num2str(M12bar(28))])
% disp(['relabel   = ',num2str(tend.fluid(1).mvareta.relabel(28)),' ',...
%                      num2str(tend.fluid(2).mvareta.relabel(28))])
% junk1 = - M21bar.*(2*Rvareta1 - 2*eta1etahat21);
% junk2 =   M21bar.*(vareta1 - vareta2 + Rvareta2 + Rvareta1 - 2*eta2etahat21);
% disp(['relabel21 = ',num2str(junk1(28)),' ',...
%                      num2str(junk2(28))])
% junk1 =   M12bar.*(vareta2 - vareta1 + Rvareta1 + Rvareta2 - 2*eta1etahat12);
% junk2 = - M12bar.*(2*Rvareta2 - 2*eta2etahat12);
% disp(['relabel12 = ',num2str(junk1(28)),' ',...
%                      num2str(junk2(28))])
% disp(['Rvareta2  = ',num2str(Rvareta2(28)),'  eta2etahat12 = ',...
%                      num2str(eta2etahat12(28))])
% Entrained and detrained values of q variance
% Note these are quasi-advective form tendencies
% Resolved variance
Rvarq1 = q1.*q1;
Rvarq2 = q2.*q2;
% 'Upwind' approximation for transferred variances
% qvarhat12 = varq1;
% qvarhat21 = varq2;
% q times transferred q
%q1qhat12 = q1.*relabel.qhat12;
%q2qhat12 = q2.*relabel.qhat12;
%q2qhat21 = q2.*relabel.qhat21;
%q1qhat21 = q1.*relabel.qhat21;
% q variance tendencies due to entrainment/detrainment                       
% tend.fluid(1).mvarq.relabel = M12bar.*(varq2 - varq1 ...
%                                      + Rvarq1 + Rvarq2 - 2*q1qhat12) ...
%                             - M21bar.*(2*Rvarq1 - 2*q1qhat21);
% tend.fluid(2).mvarq.relabel = M21bar.*(varq1 - varq2 ...
%                                      + Rvarq2 + Rvarq1 - 2*q2qhat21) ...
%                             - M12bar.*(2*Rvarq2 - 2*q2qhat12);
tend.fluid(1).mvarq.relabel = M12bar.*(varq2 - varq1 ...
                                      + (relabel.qhat12 - q1).^2) ...
                              - M21bar.*(relabel.qhat21 - q1).^2;
tend.fluid(2).mvarq.relabel = M21bar.*(varq1 - varq2 ...
                                      + (relabel.qhat21 - q2).^2) ...
                              - M12bar.*(relabel.qhat12 - q2).^2;

% ------

% Total tendencies

% Mass
tend.fluid(1).m.tot = tend.fluid(1).m.transport ...
                    + tend.fluid(1).m.relabel;
tend.fluid(2).m.tot = tend.fluid(2).m.transport ...
                    + tend.fluid(2).m.relabel;
% disp(['transport ',num2str(tend.fluid(2).m.transport(53))])
% disp(['relabel   ',num2str(tend.fluid(2).m.relabel(53))])
% disp(['M12       ',num2str(M12(53))])
% disp(['M21       ',num2str(M21(53))])


% Entropy
tend.fluid(1).meta.tot = tend.fluid(1).meta.transport ...
                       + tend.fluid(1).meta.diffuse ...
                       + tend.fluid(1).meta.diffent ...
                       + tend.fluid(1).meta.relabel ...
                       + tend.fluid(1).meta.bflux ...
                       + tend.fluid(1).meta.dissn;
tend.fluid(2).meta.tot = tend.fluid(2).meta.transport ...
                       + tend.fluid(2).meta.diffuse ...
                       + tend.fluid(2).meta.diffent ...
                       + tend.fluid(2).meta.relabel ...
                       + tend.fluid(2).meta.bflux ...
                       + tend.fluid(2).meta.dissn;

% Water
tend.fluid(1).mq.tot = tend.fluid(1).mq.transport ...
                     + tend.fluid(1).mq.diffuse ...
                     + tend.fluid(1).mq.diffent ...
                     + tend.fluid(1).mq.relabel;
tend.fluid(2).mq.tot = tend.fluid(2).mq.transport ...
                     + tend.fluid(2).mq.diffuse ...
                     + tend.fluid(2).mq.diffent ...
                     + tend.fluid(2).mq.relabel;
                 
% Vertical velocity
tend.fluid(1).mw.tot = tend.fluid(1).mw.transport ...
                     + tend.fluid(1).mw.diffuse ...
                     + tend.fluid(1).mw.diffent ...
                     + tend.fluid(1).mw.drag ...
                     + tend.fluid(1).mw.pgterm ...
                     + tend.fluid(1).mw.relabel;
tend.fluid(2).mw.tot = tend.fluid(2).mw.transport ...
                     + tend.fluid(2).mw.diffuse ...
                     + tend.fluid(2).mw.diffent ...
                     + tend.fluid(2).mw.drag ...
                     + tend.fluid(2).mw.pgterm ...
                     + tend.fluid(2).mw.relabel;
% Zero tendencies at boundaries
tend.fluid(1).mw.tot(1) = 0;
tend.fluid(1).mw.tot(nzp) = 0;
tend.fluid(2).mw.tot(1) = 0;
tend.fluid(2).mw.tot(nzp) = 0;

% Horizontal velocity
tend.fluid(1).mu.tot = tend.fluid(1).mu.transport ...
                     + tend.fluid(1).mu.diffuse ...
                     + tend.fluid(1).mu.diffent ...
                     + tend.fluid(1).mu.coriolis ...
                     + tend.fluid(1).mu.relabel;
tend.fluid(1).mv.tot = tend.fluid(1).mv.transport ...
                     + tend.fluid(1).mv.diffuse ...
                     + tend.fluid(1).mv.diffent ...
                     + tend.fluid(1).mv.coriolis ...
                     + tend.fluid(1).mv.relabel;
tend.fluid(2).mu.tot = tend.fluid(2).mu.transport ...
                     + tend.fluid(2).mu.diffuse ...
                     + tend.fluid(2).mu.diffent ...
                     + tend.fluid(2).mu.coriolis ...
                     + tend.fluid(2).mu.relabel;
tend.fluid(2).mv.tot = tend.fluid(2).mv.transport ...
                     + tend.fluid(2).mv.diffuse ...
                     + tend.fluid(2).mv.diffent ...
                     + tend.fluid(2).mv.coriolis ...
                     + tend.fluid(2).mv.relabel;

% TKE
tend.fluid(1).mtke.tot = tend.fluid(1).mtke.transport ...
                       + tend.fluid(1).mtke.diffuse ...
                       + tend.fluid(1).mtke.diffent ...
                       + tend.fluid(1).mtke.shear ...
                       + tend.fluid(1).mtke.bflux ...
                       + tend.fluid(1).mtke.drag ...
                       + tend.fluid(1).mtke.dissn ...
                       + tend.fluid(1).mtke.relabel;
tend.fluid(2).mtke.tot = tend.fluid(2).mtke.transport ...
                       + tend.fluid(2).mtke.diffuse ...
                       + tend.fluid(2).mtke.diffent ...
                       + tend.fluid(2).mtke.shear ...
                       + tend.fluid(2).mtke.bflux ...
                       + tend.fluid(2).mtke.drag ...
                       + tend.fluid(2).mtke.dissn ...
                       + tend.fluid(2).mtke.relabel;

% eta variance
tend.fluid(1).mvareta.tot = tend.fluid(1).mvareta.diffuse ...
                          + tend.fluid(1).mvareta.diffent ...
                          + tend.fluid(1).mvareta.dissn ...
                          + tend.fluid(1).mvareta.relabel;
tend.fluid(2).mvareta.tot = tend.fluid(2).mvareta.diffuse ...
                          + tend.fluid(2).mvareta.diffent ...
                          + tend.fluid(2).mvareta.dissn ...
                          + tend.fluid(2).mvareta.relabel;
                      
% q variance
tend.fluid(1).mvarq.tot = tend.fluid(1).mvarq.diffuse ...
                        + tend.fluid(1).mvarq.diffent ...
                        + tend.fluid(1).mvarq.dissn ...
                        + tend.fluid(1).mvarq.relabel;
tend.fluid(2).mvarq.tot = tend.fluid(2).mvarq.diffuse ...
                        + tend.fluid(2).mvarq.diffent ...
                        + tend.fluid(2).mvarq.dissn ...
                        + tend.fluid(2).mvarq.relabel;
                    
% For testing:
% check_equal_tendencies

% Capture advective-form budgets for diagnostics
m1transbar = weight_to_w(grid,tend.fluid(1).m.transport);
m1totbar   = weight_to_w(grid,tend.fluid(1).m.tot);
m2transbar = weight_to_w(grid,tend.fluid(2).m.transport);
m2totbar   = weight_to_w(grid,tend.fluid(2).m.tot);

budgets.w2.transport = (tend.fluid(2).mw.transport - m2transbar.*w2)./m2bar;
budgets.w2.diffuse   = (tend.fluid(2).mw.diffuse + tend.fluid(2).mw.diffent)./m2bar;
budgets.w2.entrain   =  relabel.M21bar.*(relabel.what21 - w2)./m2bar;
budgets.w2.detrain   = -relabel.M12bar.*(relabel.what12 - w2)./m2bar;
budgets.w2.pgterm    = tend.fluid(2).mw.pgterm./m2bar;
budgets.w2.drag      = tend.fluid(2).mw.drag./m2bar;
budgets.w2.tot       = (tend.fluid(2).mw.tot - m2totbar.*w2)./m2bar;

budgets.eta1.transport = (tend.fluid(1).meta.transport - m1transbar.*eta1)./m1bar;
budgets.eta1.diffuse   = (tend.fluid(1).meta.diffuse + tend.fluid(1).meta.diffent)./m1bar;
budgets.eta1.entrain   = -relabel.M21bar.*(relabel.etahat21 - eta1)./m1bar;
budgets.eta1.detrain   =  relabel.M12bar.*(relabel.etahat12 - eta1)./m1bar;
budgets.eta1.bflux     = tend.fluid(1).meta.bflux./m1bar;
budgets.eta1.dissn     = tend.fluid(1).meta.dissn./m1bar;
budgets.eta1.tot       = (tend.fluid(1).meta.tot - m1totbar.*eta1)./m1bar;

budgets.eta2.transport = (tend.fluid(2).meta.transport - m2transbar.*eta2)./m2bar;
budgets.eta2.diffuse   = (tend.fluid(2).meta.diffuse + tend.fluid(2).meta.diffent)./m2bar;
budgets.eta2.entrain   =  relabel.M21bar.*(relabel.etahat21 - eta2)./m2bar;
budgets.eta2.detrain   = -relabel.M12bar.*(relabel.etahat12 - eta2)./m2bar;
budgets.eta2.bflux     = tend.fluid(2).meta.bflux./m2bar;
budgets.eta2.dissn     = tend.fluid(2).meta.dissn./m2bar;
budgets.eta2.tot       = (tend.fluid(2).meta.tot - m2totbar.*eta2)./m2bar;

budgets.q1.transport = (tend.fluid(1).mq.transport - m1transbar.*q1)./m1bar;
budgets.q1.diffuse   = (tend.fluid(1).mq.diffuse + tend.fluid(1).mq.diffent)./m1bar;
budgets.q1.entrain   = -relabel.M21bar.*(relabel.qhat21 - q1)./m1bar;
budgets.q1.detrain   =  relabel.M12bar.*(relabel.qhat12 - q1)./m1bar;
budgets.q1.tot       = (tend.fluid(1).mq.tot - m1totbar.*q1)./m1bar;

budgets.q2.transport = (tend.fluid(2).mq.transport - m2transbar.*q2)./m2bar;
budgets.q2.diffuse   = (tend.fluid(2).mq.diffuse + tend.fluid(2).mq.diffent)./m2bar;
budgets.q2.entrain   =  relabel.M21bar.*(relabel.qhat21 - q2)./m2bar;
budgets.q2.detrain   = -relabel.M12bar.*(relabel.qhat12 - q2)./m2bar;
budgets.q2.tot       = (tend.fluid(2).mq.tot - m2totbar.*q2)./m2bar;

% Save working fields for use elsewhere
work.F1 = F1;
work.F2 = F2;
work.m1bar = m1bar;
work.m2bar = m2bar;
work.kdifft1 = kdifft1;
work.kdifft2 = kdifft2;
work.kdiffq1 = kdiffq1;
work.kdiffq2 = kdiffq2;
work.kdiffw1 = kdiffw1;
work.kdiffw2 = kdiffw2;
work.kdiffu1 = kdiffu1;
work.kdiffu2 = kdiffu2;
work.kdifftke1 = kdifftke1;
work.kdifftke2 = kdifftke2;
  work.kdifft1x = kdifft1x;
  work.kdifft2x = kdifft2x;
work.F1bar = F1bar;
work.F2bar = F2bar;
work.dF1dma = dF1dma;
work.dF1dmb = dF1dmb;
work.dF1dw = dF1dw;
work.dF2dma = dF2dma;
work.dF2dmb = dF2dmb;
work.dF2dw = dF2dw;
work.Feta1 = Feta1;
work.dFeta1detaa = dFeta1detaa;
work.dFeta1detab = dFeta1detab;
work.deta1udz = deta1udz;
work.eta1ubar = eta1ubar;
work.Feta2 = Feta2;
work.dFeta2detaa = dFeta2detaa;
work.dFeta2detab = dFeta2detab;
work.deta2udz = deta2udz;
work.eta2ubar = eta2ubar;
work.Fq1 = Fq1;
work.dFq1dqa = dFq1dqa;
work.dFq1dqb = dFq1dqb;
work.dq1udz = dq1udz;
work.q1ubar = q1ubar;
work.Fq2 = Fq2;
work.dFq2dqa = dFq2dqa;
work.dFq2dqb = dFq2dqb;
work.dq2udz = dq2udz;
work.q2ubar = q2ubar;
work.dpdz = dpdz;
work.ddragdw1 = ddragdw1;
work.ddragdw2 = ddragdw2;
work.Fw1 = Fw1;
work.dFw1dwa = dFw1dwa;
work.dFw1dwb = dFw1dwb;
work.dw1udz = dw1udz;
work.w1ubar = w1ubar;
work.Fw2 = Fw2;
work.dFw2dwa = dFw2dwa;
work.dFw2dwb = dFw2dwb;
work.dw2udz = dw2udz;
work.w2ubar = w2ubar;
work.Deta1 = Deta1;
work.Deta2 = Deta2;
work.Deta1bc = Deta1bc;
work.Deta2bc = Deta2bc;
work.Dq1 = Dq1;
work.Dq2 = Dq2;
work.Fu1 = Fu1;
work.dFu1dua = dFu1dua;
work.dFu1dub = dFu1dub;
work.du1udz = du1udz;
work.u1ubar = u1ubar;
work.Fu2 = Fu2;
work.dFu2dua = dFu2dua;
work.dFu2dub = dFu2dub;
work.du2udz = du2udz;
work.u2ubar = u2ubar;
work.Fv1 = Fv1;
work.dFv1dva = dFv1dva;
work.dFv1dvb = dFv1dvb;
work.dv1udz = dv1udz;
work.v1ubar = v1ubar;
work.Fv2 = Fv2;
work.dFv2dva = dFv2dva;
work.dFv2dvb = dFv2dvb;
work.dv2udz = dv2udz;
work.v2ubar = v2ubar;
work.dDu1dua = dDu1dua;
work.dDu1dub = dDu1dub;
work.dDv1dva = dDv1dva;
work.dDv1dvb = dDv1dvb;
work.dDu2dua = dDu2dua;
work.dDu2dub = dDu2dub;
work.dDv2dva = dDv2dva;
work.dDv2dvb = dDv2dvb;
%work.duhat12du1 = duhat12du1;
%work.duhat12du2 = duhat12du2;
%work.duhat21du1 = duhat21du1;
%work.duhat21du2 = duhat21du2;
%work.dvhat12dv1 = dvhat12dv1;
%work.dvhat12dv2 = dvhat12dv2;
%work.dvhat21dv1 = dvhat21dv1;
%work.dvhat21dv2 = dvhat21dv2;
work.dDeta1detaa = dDeta1detaa;
work.dDeta1detab = dDeta1detab;
work.dDeta2detaa = dDeta2detaa;
work.dDeta2detab = dDeta2detab;
work.dDq1dqa = dDq1dqa;
work.dDq1dqb = dDq1dqb;
work.dDq2dqa = dDq2dqa;
work.dDq2dqb = dDq2dqb;
work.dDw1dwa = dDw1dwa;
work.dDw1dwb = dDw1dwb;
work.dDw2dwa = dDw2dwa;
work.dDw2dwb = dDw2dwb;
work.dDeta1dm = dDeta1dm;
work.dDeta2dm = dDeta2dm;
work.dDq1dm = dDq1dm;
work.dDq2dm = dDq2dm;
work.dDw1dm = dDw1dm;
work.dDw2dm = dDw2dm;
work.nhpg1 = nhpg1;
work.nhpg2 = nhpg2;
% work.L_turb1 = scales.L_turb1;
% work.L_turb2 = scales.L_turb2;
work.dFtke1dtkea = dFtke1dtkea;
work.dFtke1dtkeb = dFtke1dtkeb;
work.dtke1udz = dtke1udz;
work.tke1ubar = tke1ubar;
work.dFtke2dtkea = dFtke2dtkea;
work.dFtke2dtkeb = dFtke2dtkeb;
work.dtke2udz = dtke2udz;
work.tke2ubar = tke2ubar;
work.dDtke1dtkea = dDtke1dtkea;
work.dDtke1dtkeb = dDtke1dtkeb;
work.dDtke2dtkea = dDtke2dtkea;
work.dDtke2dtkeb = dDtke2dtkeb;
%work.dLdtke1 = dLdtke1;
%work.dLdtke2 = dLdtke2;
work.buoy = buoy;
%work.buoyx = buoyx;
work.bflux1 = bflux1;
work.bflux2 = bflux2;
work.bflux_res = buoy_flux_res;


% ------

end