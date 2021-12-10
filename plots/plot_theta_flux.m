% Plot contributions to theta flux

% In the compressible case there can be a non-trivial contribution
% from mean ascent due to warming and expansion of the boundary layer
% air; we should be careful to separate this out.

% ----------

% Figure 29: contributions to theta flux

figure(29)

% Mean density on w-levels
rhow = m1bar + m2bar;

% Mass weighted mean w
w_bar_star = (work.F1 + work.F2)./rhow;

% Mass weighted mean theta
theta_bar_star = (m1bar.*eos.theta1 + m2bar.*eos.theta2)./rhow;

% Deviation mass fluxes
F1_pert = work.F1 - m1bar.*w_bar_star;
F2_pert = work.F2 - m2bar.*w_bar_star;

% Mean ascent contribution
HF_mean = rhow.*w_bar_star.*theta_bar_star;

% Resolved contributions in fluids 1 and 2
HF_res1 = F1_pert.*(eos.theta1 - theta_bar_star);
HF_res2 = F2_pert.*(eos.theta2 - theta_bar_star);

% Subfilter-scale contributions
HF_sf1 = work.bflux1.*eos.theta1p/settings.constants.phys.gravity;
HF_sf2 = work.bflux2.*eos.theta2p/settings.constants.phys.gravity;

% Total
HF_tot = HF_mean + HF_res1 + HF_res2 + weight_to_w(grid,HF_sf1 + HF_sf2);

% And plot

subplot(1,3,1)
plot(HF_tot,zunitsw,'k',...
     HF_mean,zunitsw,'r')
ylim([0,zplottop])
title('Tot and expansion')
xlabel('w \theta')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(1,3,2)
plot(HF_res1,zunitsw,'b',...
     HF_res2,zunitsw,'r')
ylim([0,zplottop])
title('Resolved')
xlabel('w \theta')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(1,3,3)
plot(HF_sf1,zunitsp,'b',...
     HF_sf2,zunitsp,'r')
ylim([0,zplottop])
title('Subfilter')
xlabel('w \theta')
ylabel(labelz)
set(gca,'FontSize',fs)

