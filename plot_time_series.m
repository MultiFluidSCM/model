% Plot snap shot profiles and various time series

% Figure 2

% List of times at which to plot theta profile
profile_times = [0, 10800, 21600, 32400, 43200];
% profile_times = [0, 14000, 19400, 28400, 33200];
npt = numel(profile_times); 

% Time step size
dt = time.dt;


% Plot theta profile at selected times
% Is it time to plot?
if plotstep == 1
    lasttime = -1;
else
    lasttime = ts.time(plotstep-1);
end
lprofile = 0;
ipt = 0;
while ipt < npt & ~lprofile
    ipt = ipt + 1;
    pt = profile_times(ipt);
    if (lasttime < pt & time.t >= pt)
        lprofile = 1;
    end
end
if lprofile
    % Mean potential temperature
    % thetaave = (m1bar.*eos.theta_rho1 + m2bar.*eos.theta_rho2)./(m1bar + m2bar);
    thetaave = (m1bar.*eos.theta1 + m2bar.*eos.theta2)./(m1bar + m2bar);
    if plottype == 0
        figure(2)
        subplot(2,3,1)
        if ipt == 1
            hold off
        end
        plot(thetaave,grid.zwkm,'k')
        axis([297,315,0,4])
        title('Theta')
        xlabel('\theta_l (K)')
        ylabel('Height (km)')
        set(gca,'FontSize',fs)
        hold on
    % else
        figure(8)
        subplot(1,1,1)
        if ipt == 1
            hold off
        end
        plot(thetaave,grid.zwkm,'k','LineWidth',1.0)
        % axis([297,315,0,3])
        axis([295,320,0,4])
        title('Theta')
        xlabel('\theta_l (K)')
        ylabel('Height (km)')
        set(gca,'FontSize',fs)
        hold on
    end
end


% Time series of zstar
ts.zstar(plotstep) = zstar;
if plottype == 0
    figure(2)
    subplot(2,3,4)
    plot(ts.time/3600,ts.zstar,'k')
    axis([0,time.tstop/3600,0,grid.zw(nzp)])
    title(['t = ',num2str(time.t),'  z* = ',num2str(zstar,'%6.1f')])
    xlabel('time')
    ylabel('z_*, z_{cb}, z_{ct}')
    set(gca,'FontSize',fs)
else
    figure(6)
    subplot(3,1,2)
    plot(ts.time/3600,ts.zstar,'k','linewidth',1.5)
    xlim([0,52240/3600])
    title(['t = ',num2str(time.t)])
    xlabel('time')
    ylabel('z_*, z_{cb}, z_{ct}')
    set(gca,'FontSize',fs)
end
% and time series of cloud base and top
zcldbase = 0;
k = 1;
lq = liquid1(k) + liquid2(k);
while lq < 1e-5 & k < nzp
    k = k + 1;
    lq = liquid1(k) + liquid2(k);
    if lq > 1e-5
        zcldbase = grid.zw(k);
    end
end
ts.zcbase(plotstep) = zcldbase;
ts.zcbaseSG(plotstep) = z_cld_base;
zcldtop = 0;
k = nzp;
lq = liquid1(k) + liquid2(k);
while lq < 1e-5 & k > 1
    k = k - 1;
    lq = liquid1(k) + liquid2(k);
    if lq > 1e-5
        zcldtop = grid.zw(k);
    end
end
ts.zctop(plotstep) = zcldtop;
ts.zctopSG(plotstep) = z_cld_top;
ts.lnbgas(plotstep) = scales.LNBgas;
if plottype == 0
    figure(2)
    subplot(2,3,4)
    hold on
    %plot(ts.time/3600,ts.zcbase,'r--')
    %plot(ts.time/3600,ts.zctop, 'b--')
    plot(ts.time/3600,ts.zcbaseSG,'r')
    plot(ts.time/3600,ts.zctopSG ,'b')
    % plot(ts.time/3600,ts.lnbgas,'g')
    hold off
else
    figure(6)
    subplot(3,1,2)
    hold on
    plot(ts.time/3600,ts.zcbaseSG,'r','linewidth',1.5)
    plot(ts.time/3600,ts.zctopSG ,'b','linewidth',1.5)
    % plot(ts.time/3600,ts.lnbgas,'g','linewidth',1.5)
    hold off
    pause
end
%ts.tke_thresh(plotstep) = scales.tke_thresh;
ts.totcldcov(plotstep) = tot_cld_cov;
figure(21)
subplot(1,3,3)
plot(ts.time/3600,ts.totcldcov,'k','linewidth',1.2)
%xlim([0,time.tstop/3600])
axis([0,time.tstop/3600,0,0.4])
set(gca,'FontSize',fs)
title(['t = ',num2str(time.t),'  z* = ',num2str(zstar,'%6.1f')])
xlabel('time')
ylabel('Cloud cover')


% % Time series of column liquid water
% colliquid = sum(m1bar.*liquid1 + m2bar.*liquid2);
% ts.liquid(plotstep) = colliquid;
% if plottype > 0
%     figure(4)
%     subplot(3,1,3)
%     plot(ts.time/3600,ts.liquid,'k','linewidth',1.5)
%     %xlim([0,time.tstop/3600])
%     xlim([0,52240/3600])
%     title(['t = ',num2str(time.t)])
%     xlabel('time')
%     ylabel('Col. liquid (kg/m^2)')
%     set(gca,'FontSize',fs)
% end


% Time series of total mass, entropy, water, and
% and energy change
ts.mass(plotstep) = gdiags.mass1 + gdiags.mass2;
mass0 = ts.mass(1);
ts.entropy(plotstep) = gdiags.entropy1 + gdiags.entropy2;
entropy0 = ts.entropy(1);
ts.water(plotstep) = gdiags.water1 + gdiags.water2;
water0 = ts.water(1);
ts.energy(plotstep) = gdiags.energy1 + gdiags.energy2;
energy0 = ts.energy(1);
ts.mass_force(plotstep)    = accum_force.smf;
ts.entropy_force(plotstep) = accum_force.setaf;
ts.water_force(plotstep)   = accum_force.sqf;
ts.energy_force(plotstep)  = accum_force.sEf;
ts.etabflux(plotstep)      = accum.eta_bflux;
ts.etadissn(plotstep)      = accum.eta_dissn;
ts.dtke_fix(plotstep)      = accum.tke1_fix + accum.tke2_fix;

if plottype == 0
    figure(2)
    subplot(2,6,3)
    plot(ts.time/3600,(ts.mass - mass0)/mass0,'b',...
         ts.time/3600,(ts.entropy - entropy0)/entropy0,'r',...
         ts.time/3600,(ts.water - water0)/water0,'g',...
         ts.time/3600,(ts.energy - energy0)/energy0,'k')
    xlim([0,time.tstop/3600])
    title('Rel. chg.')
    xlabel('time')
    ylabel('Rel. chg.')
    legend('M','S','Q','E','Location','NorthWest')
    set(gca,'FontSize',fs)
    subplot(2,6,9)
    plot(ts.time/3600,(ts.mass - mass0 - ts.mass_force)/mass0,'b',...
         ts.time/3600,(ts.entropy - entropy0 - ts.entropy_force - ts.etabflux - ts.etadissn)/entropy0,'r',...
         ts.time/3600,(ts.water - water0 - ts.water_force)/water0,'g',...
         ts.time/3600,(ts.energy - energy0 - ts.energy_force - ts.dtke_fix)/energy0,'k')
    xlim([0,time.tstop/3600])
    title('Discrepancy')
    xlabel('time')
    ylabel('Discrepancy')
    set(gca,'FontSize',fs)
else
    figure(7)
    subplot(2,4,1)
    plot(ts.time/3600,(ts.mass - mass0)/mass0,'r')
    hold on
    plot(ts.time/3600,ts.mass_force/mass0,'b')
    hold off
    xlim([0,time.tstop/3600])
    xlabel('time')
    ylabel('Rel. chg.')
    title('Mass change')
    set(gca,'FontSize',fs)
    subplot(2,4,2)
    plot(ts.time/3600,(ts.entropy - entropy0)/entropy0,'r')
    hold on
    plot(ts.time/3600,ts.entropy_force/entropy0,'b')
    hold off
    xlim([0,time.tstop/3600])
    xlabel('time')
    ylabel('Rel. chg.')
    title('Entropy change')
    set(gca,'FontSize',fs)
    subplot(2,4,3)
    plot(ts.time/3600,(ts.water - water0)/water0,'r')
    hold on
    plot(ts.time/3600,ts.water_force/water0,'b')
    hold off
    xlim([0,time.tstop/3600])
    xlabel('time')
    ylabel('Rel. chg.')
    title('Water change')
    set(gca,'FontSize',fs)
    subplot(2,4,4)
    set(gca,'FontSize',fs)
    plot(ts.time/3600,(ts.energy - energy0)/energy0,'r')
    hold on
    plot(ts.time/3600,ts.energy_force/energy0,'b')
    hold off
    xlim([0,time.tstop/3600])
    xlabel('time')
    ylabel('Rel. chg.')
    title('Energy change')
    subplot(2,4,5)
    plot(ts.time/3600,(ts.mass - mass0 - ts.mass_force)/mass0,'r')
    xlim([0,time.tstop/3600])
    xlabel('time')
    ylabel('Disc.')
    title('Mass discrepancy')
    set(gca,'FontSize',fs)
    subplot(2,4,6)
    plot(ts.time/3600,(ts.entropy - entropy0 - ts.entropy_force - ts.etabflux - ts.etadissn)/entropy0,'r')
    xlim([0,time.tstop/3600])
    xlabel('time')
    ylabel('Disc.')
    title('Entropy discrepancy')
    set(gca,'FontSize',fs)
    subplot(2,4,7)
    plot(ts.time/3600,(ts.water - water0 - ts.water_force)/water0,'r')
    xlim([0,time.tstop/3600])
    xlabel('time')
    ylabel('Disc.')
    title('Water discrepancy')
    set(gca,'FontSize',fs)
    subplot(2,4,8)
    plot(ts.time/3600,(ts.energy - energy0 - ts.energy_force - ts.dtke_fix)/energy0,'r')
    xlim([0,time.tstop/3600])
    xlabel('time')
    ylabel('Disc.')
    title('Energy discrepancy')
    set(gca,'FontSize',fs)
    pause
end


% Accumulated KE budgets
ts.dke1_pg(plotstep)      = accum.ke1_pg;
ts.dke1_coriol(plotstep)  = accum.ke1_coriol;
ts.dke1_diff(plotstep)    = accum.ke1_diff;
ts.dke1_drag(plotstep)    = accum.ke1_drag;
ts.dke1_relabel(plotstep) = accum.ke1_relabel;
ts.dke2_pg(plotstep)      = accum.ke2_pg;
ts.dke2_coriol(plotstep)  = accum.ke2_coriol;
ts.dke2_diff(plotstep)    = accum.ke2_diff;
ts.dke2_drag(plotstep)    = accum.ke2_drag;
ts.dke2_relabel(plotstep) = accum.ke2_relabel;
ts.dke1_tot = ts.dke1_pg + ts.dke1_coriol + ts.dke1_diff + ts.dke1_drag + ts.dke1_relabel;
ts.dke2_tot = ts.dke2_pg + ts.dke2_coriol + ts.dke2_diff + ts.dke2_drag + ts.dke2_relabel;
ts.ke1_check(plotstep) = gdiags.kinetic1;
ts.ke2_check(plotstep) = gdiags.kinetic2;
if plottype == 0
    figure(2)
    subplot(2,6,4)
    plot(ts.time/3600,ts.dke1_pg,'b',...
         ts.time/3600,ts.dke1_coriol,'r--',...
         ts.time/3600,ts.dke1_diff,'r',...
         ts.time/3600,ts.dke1_drag,'g',...
         ts.time/3600,ts.dke1_relabel,'b--',...
         ts.time/3600,ts.dke1_tot,'k',...
         ts.time/3600,ts.ke1_check,'k')
    xlim([0,time.tstop/3600])
    title('Cum. chng. in res KE1')
    xlabel('time')
    %ylabel('Dke1 PG, C, K, P, R')
    legend('PG','C','K','P','KH','R','T','Location','SouthWest')
    set(gca,'FontSize',fs)
    subplot(2,6,5)
    plot(ts.time/3600,ts.dke2_pg,'b',...
         ts.time/3600,ts.dke2_coriol,'r--',...
         ts.time/3600,ts.dke2_diff,'r',...
         ts.time/3600,ts.dke2_drag,'g',...
         ts.time/3600,ts.dke2_relabel,'b--',...
         ts.time/3600,ts.dke2_tot,'k',...
         ts.time/3600,ts.ke2_check,'k')
    xlim([0,time.tstop/3600])
    title('Cum. chng. in res KE2')
    xlabel('time')
    %ylabel('Dke2 PG, C, K, P, R')
    set(gca,'FontSize',fs)
    subplot(2,6,6)
    plot(ts.time/3600,ts.dke1_pg     +ts.dke2_pg,'b',...
         ts.time/3600,ts.dke1_coriol +ts.dke2_coriol,'r--',...
         ts.time/3600,ts.dke1_diff   +ts.dke2_diff,'r',...
         ts.time/3600,ts.dke1_drag   +ts.dke2_drag,'g',...
         ts.time/3600,ts.dke1_relabel+ts.dke2_relabel,'b--',...
         ts.time/3600,ts.dke1_tot    +ts.dke2_tot,'k')
    xlim([0,time.tstop/3600])
    title('Cum. chng. in res KE')
    xlabel('time')
    %ylabel('Dke PG, C, K, P, R')
    set(gca,'FontSize',fs)
else
    figure(7)
    subplot(2,3,1)
    plot(ts.time/3600,ts.dke1_pg,'b',...
         ts.time/3600,ts.dke1_coriol,'r--',...
         ts.time/3600,ts.dke1_diff,'r',...
         ts.time/3600,ts.dke1_drag,'g',...
         ts.time/3600,ts.dke1_relabel,'b--',...
         ts.time/3600,ts.dke1_tot,'k')
    xlim([0,time.tstop/3600])
    title('Cum. chng. in res KE1')
    xlabel('time')
    %ylabel('Dke1 PG, C, K, P, R')
    legend('PG','C','K','P','R','T','Location','SouthWest')
    set(gca,'FontSize',fs)
    subplot(2,3,2)
    plot(ts.time/3600,ts.dke2_pg,'b',...
         ts.time/3600,ts.dke2_coriol,'r--',...
         ts.time/3600,ts.dke2_diff,'r',...
         ts.time/3600,ts.dke2_drag,'g',...
         ts.time/3600,ts.dke2_relabel,'b--',...
         ts.time/3600,ts.dke2_tot,'k')
    xlim([0,time.tstop/3600])
    title('Cum. chng. in res KE2')
    xlabel('time')
    %ylabel('Dke2 PG, C, K, P, R')
    set(gca,'FontSize',fs)
    subplot(2,3,3)
    plot(ts.time/3600,ts.dke1_pg     +ts.dke2_pg,'b',...
         ts.time/3600,ts.dke1_coriol +ts.dke2_coriol,'r--',...
         ts.time/3600,ts.dke1_diff   +ts.dke2_diff,'r',...
         ts.time/3600,ts.dke1_drag   +ts.dke2_drag,'g',...
         ts.time/3600,ts.dke1_relabel+ts.dke2_relabel,'b--',...
         ts.time/3600,ts.dke1_tot    +ts.dke2_tot,'k')
    xlim([0,time.tstop/3600])
    title('Cum. chng. in res KE')
    xlabel('time')
    %ylabel('Dke PG, C, K, P, R')
    set(gca,'FontSize',fs)
end
    

% Accumulated TKE budgets
ts.tke1(plotstep) = gdiags.tke1;
ts.tke2(plotstep) = gdiags.tke2;
tke0 = ts.tke1(1) + ts.tke2(1);
ts.dtke1_shear(plotstep)   = accum.tke1_shear;
ts.dtke1_buoy(plotstep)    = accum.tke1_buoy;
ts.dtke1_drag(plotstep)    = accum.tke1_drag;
ts.dtke1_diss(plotstep)    = accum.tke1_diss;
ts.dtke1_fix(plotstep)     = accum.tke1_fix;
ts.dtke1_relabel(plotstep) = accum.tke1_relabel;
ts.dtke2_shear(plotstep)   = accum.tke2_shear;
ts.dtke2_buoy(plotstep)    = accum.tke2_buoy;
ts.dtke2_drag(plotstep)    = accum.tke2_drag;
ts.dtke2_diss(plotstep)    = accum.tke2_diss;
ts.dtke2_fix(plotstep)     = accum.tke2_fix;
ts.dtke2_relabel(plotstep) = accum.tke2_relabel;
ts.dtke1_tot = ts.dtke1_shear + ts.dtke1_buoy    + ts.dtke1_drag + ...
               ts.dtke1_diss  + ts.dtke1_relabel + ts.dtke1_fix;
ts.dtke2_tot = ts.dtke2_shear + ts.dtke2_buoy    + ts.dtke2_drag ...
             + ts.dtke2_diss  + ts.dtke2_relabel + ts.dtke2_fix;
if plottype == 0
    figure(2)
    subplot(2,6,10)
    plot(ts.time/3600,ts.dtke1_shear,'b',...
         ts.time/3600,ts.dtke1_buoy,'r--',...
         ts.time/3600,ts.dtke1_drag,'r',...
         ts.time/3600,ts.dtke1_diss,'g',...
         ts.time/3600,ts.dtke1_relabel,'b--',...
         ts.time/3600,ts.dtke1_fix,'k:',...
         ts.time/3600,ts.dtke1_tot,'k',...
         ts.time/3600,ts.tke1,'k--')
    xlim([0,time.tstop/3600])
    title('Cum. chng. in TKE1')
    xlabel('time')
    legend('sh','B','P','D','R','X','T','A','Location','SouthWest')
    set(gca,'FontSize',fs)
    subplot(2,6,11)
    plot(ts.time/3600,ts.dtke2_shear,'b',...
         ts.time/3600,ts.dtke2_buoy,'r--',...
         ts.time/3600,ts.dtke2_drag,'r',...
         ts.time/3600,ts.dtke2_diss,'g',...
         ts.time/3600,ts.dtke2_relabel,'b--',...
         ts.time/3600,ts.dtke2_fix,'k:',...
         ts.time/3600,ts.dtke2_tot,'k',...
         ts.time/3600,ts.tke2,'k--')
    xlim([0,time.tstop/3600])
    title('Cum. chng. in TKE2')
    xlabel('time')
    set(gca,'FontSize',fs)
    subplot(2,6,12)
    plot(ts.time/3600,ts.dtke1_shear  +ts.dtke2_shear,'b',...
         ts.time/3600,ts.dtke1_buoy   +ts.dtke2_buoy,'r--',...
         ts.time/3600,ts.dtke1_drag   +ts.dtke2_drag,'r',...
         ts.time/3600,ts.dtke1_diss   +ts.dtke2_diss,'g',...
         ts.time/3600,ts.dtke1_relabel+ts.dtke2_relabel,'b--',...
         ts.time/3600,ts.dtke1_fix    +ts.dtke2_fix,'k:',...
         ts.time/3600,ts.dtke1_tot    +ts.dtke2_tot,'k',...
         ts.time/3600,ts.tke1 + ts.tke2,'k--')
    xlim([0,time.tstop/3600])
    title('Cum. chng. in TKE')
    xlabel('time')
    %ylabel('Dke PG, C, K, P, R')
    set(gca,'FontSize',fs)
else
    figure(7)
    subplot(2,3,6)
    plot(ts.time/3600,ts.dtke1_shear,'b',...
         ts.time/3600,ts.dtke1_buoy,'r--',...
         ts.time/3600,ts.dtke1_drag,'r',...
         ts.time/3600,ts.dtke1_diss,'g',...
         ts.time/3600,ts.dtke1_relabel,'b--',...
         ts.time/3600,ts.dtke1_tot,'k',...
         ts.time/3600,ts.tke1,'k--')
    xlim([0,time.tstop/3600])
    title('Cum. chng. in TKE1')
    xlabel('time')
    legend('sh','B','P','D','R','T','A','Location','SouthWest')
    set(gca,'FontSize',fs)
    subplot(2,3,5)
    plot(ts.time/3600,ts.dtke2_shear,'b',...
         ts.time/3600,ts.dtke2_buoy,'r--',...
         ts.time/3600,ts.dtke2_drag,'r',...
         ts.time/3600,ts.dtke2_diss,'g',...
         ts.time/3600,ts.dtke2_relabel,'b--',...
         ts.time/3600,ts.dtke2_tot,'k',...
         ts.time/3600,ts.tke2,'k--')
    xlim([0,time.tstop/3600])
    title('Cum. chng. in TKE2')
    xlabel('time')
    set(gca,'FontSize',fs)
    subplot(2,3,6)
    plot(ts.time/3600,ts.dtke1_shear  +ts.dtke2_shear,'b',...
         ts.time/3600,ts.dtke1_buoy   +ts.dtke2_buoy,'r--',...
         ts.time/3600,ts.dtke1_drag   +ts.dtke2_drag,'r',...
         ts.time/3600,ts.dtke1_diss   +ts.dtke2_diss,'g',...
         ts.time/3600,ts.dtke1_relabel+ts.dtke2_relabel,'b--',...
         ts.time/3600,ts.dtke1_tot    +ts.dtke2_tot,'k',...
         ts.time/3600,ts.tke1 + ts.tke2,'k--')
    xlim([0,time.tstop/3600])
    title('Cum. chng. in TKE')
    xlabel('time')
    %ylabel('Dke PG, C, K, P, R')
    set(gca,'FontSize',fs)
end
    
if plottype == 1
    figure(7)
    subplot(2,3,4)
    plot(ts.time/3600,ts.tke1 + ts.tke2 - tke0 - (ts.dtke1_tot + ts.dtke2_tot),'r')
    xlim([0,time.tstop/3600])
    xlabel('time')
    ylabel('Disc.')
    title('TKE discrepancy')
    set(gca,'FontSize',fs)
end
    
    