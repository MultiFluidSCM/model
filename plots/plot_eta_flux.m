% Plot contributions to eta flux

% In the compressible case there can be a non-trivial contribution
% from mean ascent due to warming and expansion of the boundary layer
% air; we should be careful to separate this out.

% ----------

% Figure 32: contributions to eta flux

figure(32)

% Mean density on w-levels
rhow = m1bar + m2bar;

% Mass weighted mean w
w_bar_star = (work.F1 + work.F2)./rhow;

% Mass weighted mean eta
eta_bar_star = (m1bar.*eta1 + m2bar.*eta2)./rhow;

% Deviation mass fluxes
F1_pert = work.F1 - m1bar.*w_bar_star;
F2_pert = work.F2 - m2bar.*w_bar_star;

% Mean ascent contribution
HF_mean = rhow.*w_bar_star.*eta_bar_star;

% Resolved contributions in fluids 1 and 2
HF_res1 = F1_pert.*(eta1 - eta_bar_star);
HF_res2 = F2_pert.*(eta2 - eta_bar_star);

% Subfilter-scale contributions
HF_sf1 = Deta1;
HF_sf2 = Deta2;
HF_sf1_bc = work.Deta1bc;
HF_sf2_bc = work.Deta2bc;
HF_sf1_ed = work.Deta1ed;
HF_sf2_ed = work.Deta2ed;

% Total
HF_tot = HF_mean + HF_res1 + HF_res2 + weight_to_w(grid,HF_sf1 + HF_sf2);

% And plot

subplot(1,3,1)
plot(HF_tot,zunitsw,'k',...
     HF_mean,zunitsw,'r')
ylim([0,zplottop])
title('Tot and expansion')
xlabel('w \eta')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(1,3,2)
plot(HF_res1,zunitsw,'b',...
     HF_res2,zunitsw,'r')
ylim([0,zplottop])
title('Resolved')
xlabel('w \eta')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(1,3,3)
plot(HF_sf1_ed,zunitsp,'b',...
     HF_sf2_ed,zunitsp,'r',...
     HF_sf1_bc,zunitsp,'b--',...
     HF_sf2_bc,zunitsp,'r--')
ylim([0,zplottop])
title('Subfilter')
xlabel('w \eta')
ylabel(labelz)
set(gca,'FontSize',fs)

