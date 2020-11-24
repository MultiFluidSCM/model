% Plot a range of diagnostics

% Font size
fs = 16;

% Units for z axis on most plots
% km
zunitsp = grid.zpkm;
zunitsw = grid.zwkm;
zplottop = 4.4;
labelz = 'z (km)';
% Or non-dimensionalized
% zunitsp = grid.zp/scales.zstar;
% zunitsw = grid.zw/scales.zstar;
% zplottop = 1.5;
% labelz = 'z / z_*';


% Unpack fields
nz = grid.nz;
nzp = nz + 1;
dzp = grid.dzp;
dzw = grid.dzw;
abovep = grid.abovep;
belowp = grid.belowp;
abovew = grid.abovew;
beloww = grid.beloww;
abover = grid.abover;
belowr = grid.belowr;
m1   = state_new.fluid(1).m;
m2   = state_new.fluid(2).m;
sigma1 = m1./eos.rho1;
sigma2 = m2./eos.rho2;
w1   = state_new.fluid(1).w;
w2   = state_new.fluid(2).w;
eta1 = state_new.fluid(1).eta;
eta2 = state_new.fluid(2).eta;
q1   = state_new.fluid(1).q;
q2   = state_new.fluid(2).q;
p    = state_new.p;
T1   = state_new.fluid(1).T;
T2   = state_new.fluid(2).T;
Tw1  = state_new.fluid(1).Tw;
Tw2  = state_new.fluid(2).Tw;
u1   = state_new.fluid(1).u;
u2   = state_new.fluid(2).u;
v1   = state_new.fluid(1).v;
v2   = state_new.fluid(2).v;
tke1 = state_new.fluid(1).tke;
tke2 = state_new.fluid(2).tke;
m1bar = work.m1bar;
m2bar = work.m2bar;
F1bar = work.F1bar;
F2bar = work.F2bar;
Feta1 = work.Feta1;
Feta2 = work.Feta2;
Deta1 = work.Deta1;
Deta2 = work.Deta2;
Fq1 = work.Fq1;
Fq2 = work.Fq2;
Dq1 = work.Dq1;
Dq2 = work.Dq2;
mbar = m1bar + m2bar;
zstar = scales.zstar;
wstar = scales.wstar;
etastar = scales.etastar;
thetastar = scales.thetastar;
force2 = set_forcing(time.t);
qstar = scales.qstar;
setaf = surface_flux.eta;
sqf   = force2.sqf;

% Average to p levels
eta1bar = abovep.*eta1(2:nzp)...
        + belowp.*eta1(1:nz );
q1bar   = abovep.*q1(2:nzp)...
        + belowp.*q1(1:nz );
w1bar   = abovep.*w1(2:nzp)...
        + belowp.*w1(1:nz );
eta2bar = abovep.*eta2(2:nzp)...
        + belowp.*eta2(1:nz );
q2bar   = abovep.*q2(2:nzp)...
        + belowp.*q2(1:nz );
w2bar   = abovep.*w2(2:nzp)...
        + belowp.*w2(1:nz );
    
% Average to w levels
sigma1w = weight_to_w(grid,sigma1);
sigma2w = weight_to_w(grid,sigma2);

% ----------

% Diagnose cloud fractions
diagnose_cloud_frac

% ----------

% Figure 1: Basic fields mass frac, w, eta, q, u and v, buoyancy, N^2,
% thetal, RH, tke

plot_basic_fields

% ----------

% Figure 2: Profile snapshots and time series diagnostics

% Append current time to array for plotting time series
if exist('ts')
    plotstep = numel(ts.time) + 1;
else
    plotstep = 1;
end
ts.time(plotstep) = time.t;

plot_time_series

% ----------

% Figure 3: Budget diagnostics

% plot_budgets

% ----------

% Figure 4: Turbulence related diagnostics

% plot_turbulence

% ----------

% Figure 18: For development: actual and candidate entrainment and
% detrainment profiles

% plotED

% ----------

% Figure 19: Plot buoyancy of a selection of parcels lifted adiabatically

% plot_adiabats

% ----------

% Figure 20: Plot estimated standard deviations of updraft tracers

plot_std

% ----------

% Figure 21: Plot estimated cloud fractions using APDF

plot_cloud_frac

% ----------

% Figure 22: Compare N and turbulence inverse time scale

% plot_T_turb

% ----------

% Figure 23: Plot variance tendencies for updraft tracers

% plot_var_tend

% ----------

% Figure 24: Plot diffusive buoyancy fluxes

% plot_bflux

% ----------% 

% Figure 25: Plot diffusive and buoyancy correlation eta fluxes

% plot_bc_flux

% ----------% 

% Figure 26: Check that eta and q variance tendencies balance

% plot_var_budgets

% ----------% 

% Figure 27: Plot time series of max residuals

plot_res_time_series

% ----------% 

% Figure 28: 

plot_ql_var

% ----------% 

% 
% % Check nonhydrostatic PG balances w variance gradient
% % Crude normalization ignoring density
% wnorm = wstar*wstar/zstar;
% % calculate average w-level density
% sigma1 = m1./eos.rho1;
% sigma2 = m2./eos.rho2;
% sigma1bar(2:nz) = abovew(2:nz).*sigma1(2:nz) + beloww(2:nz).*sigma1(1:nz-1);
% sigma2bar(2:nz) = abovew(2:nz).*sigma2(2:nz) + beloww(2:nz).*sigma2(1:nz-1);
% rhow(2:nz) = sigma1bar(2:nz).*eos.rhow1(2:nz) + sigma2bar(2:nz).*eos.rhow2(2:nz);
% rhow(1) = 0;
% rhow(nzp) = 0;
% % Deviation of pressure from hydrostatic
% pgpert = work.dpdz + rhow*constants.phys.gravity;
% pgpert = pgpert/wnorm;
% % Resolved w variance
% wvar = m1bar.*w1.*w1 + m2bar.*w2.*w2;
% dwvar = (wvar(2:nzp) - wvar(1:nz))./dzp;
% dwvar = dwvar/wnorm;
% if plottype == 0
%     subplot(3,6,11)
%     set(gca,'FontSize',fs);
%     plot(dwvar,zndp,'b',pgpert,zndw,'r')
%     axis([-1.5,1.5,0,zndtop])
%     %%axis([0,0.2,0.85,1.15])
%     xlabel('wvar_z and pg-pert')
%     ylabel('z / z_*')
% else
%     figure(2)
%     set(gca,'FontSize',fs);
%     plot(dwvar,zndp,'b',pgpert,zndw,'r')
%     axis([-1.5,1.5,0,zndtop])
%     xlabel('wvar_z and pg-pert')
%     ylabel('z / z_*')
%     pause
%     figure(1)
% end
% 



% 
% % Check various flavours of Courant number
% dt = time.dt;
% cn = dt*w2(1:nzp)./dzw(1:nzp);
% kcn = dt*work.kdifft1./(dzp.*dzp);
% entcn = entrate*dt;
% detcn = detrate*dt;
% if plottype == 0
%     subplot(3,6,18)
%     set(gca,'FontSize',fs)
%     plot(cn,grid.zw,'r',kcn,grid.zp,'g',entcn,grid.zp,'r:',detcn,grid.zp,'b:')
%     ylim([0,grid.zw(nzp)])
%     xlim([0,4])
%     title('cn, kcn, ecn, dcn')
% else
%     figure(2)
%     set(gca,'FontSize',fs)
%     plot(cn,grid.zw,'r',kcn,grid.zp,'g',entcn,grid.zp,'r:',detcn,grid.zp,'b:')
%     ylim([0,grid.zw(nzp)])
%     title('cn, kcn, ecn, dcn')
%     pause
%     figure(1)
% end
% 
% 


% 
% % Time series of minimum theta flux
% % and zstar - z(k)
% tsmintf(istep+1) = -wtmin/qstar;
% tsdzstar(istep+1) = (zstar - zw(k_bltop-1))/(zw(2) - zw(1));
% % And smoothed version
% nshhour = ceil(1800/dt);
% if (istep >= 2*nshhour)
%     tssmoothtf(istep+1-nshhour) = sum(tsmintf(istep+1-2*nshhour:istep+1))/(1+2*nshhour);
% end
% if lpanel
%     subplot(3,5,13)
%     set(gca,'FontSize',fs)
%     %plot(tstime,tsmintf,'b',tstime,tsdzstar,'r')
%     %axis([0,nstop*dt,0,1])
%     plot(tstime/3600,tsmintf,'b')
%     hold on
%     if (istep >= 2*nshhour)
%         plot(tstime(nshhour:istep+1-nshhour)/3600,tssmoothtf(nshhour:istep+1-nshhour),'b:')
%     end
%     axis([0,nstop*dt/3600,0,0.3])
%     xlabel('time')
%     ylabel('Entrainment flux')
%     hold off
% 
%     subplot(3,5,14)
%     set(gca,'FontSize',fs)
%     %plot(tstime,tsmintf,'b',tstime,tsdzstar,'r')
%     %axis([0,nstop*dt,0,1])
%     semilogx(tstime*Nbv,tsmintf,'b')
%     axis([10,1000,0.05,0.2])
%     xlabel('tN')
%     ylabel('Entrainment flux')
% end
% if lindiv
%     figure(2)
%     set(gca,'FontSize',fs)
%     %plot(tstime,tsmintf,'b',tstime,tsdzstar,'r')
%     %axis([0,nstop*dt,0,1])
%     plot(tstime/3600,tsmintf,'b')
%     hold on
%     if (istep >= 2*nshhour)
%         plot(tstime(nshhour:istep+1-nshhour)/3600,tssmoothtf(nshhour:istep+1-nshhour),'b:')
%     end
%     axis([0,nstop*dt/3600,0,0.3])
%     xlabel('time')
%     ylabel('Entrainment flux')
%     hold off
%     pause
%     figure(1)
% end
% 

pause(0.01)