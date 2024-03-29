% Plot contributions to q flux

% In the compressible case there can be a non-trivial contribution
% from mean ascent due to warming and expansion of the boundary layer
% air; we should be careful to separate this out.

% ----------

% Figure 33: contributions to eta flux

figure(33)

% Mean density on w-levels
rhow = m1bar + m2bar;

% Mass weighted mean w
w_bar_star = (work.F1 + work.F2)./rhow;

% Mass weighted mean q
q_bar_star = (m1bar.*q1 + m2bar.*q2)./rhow;

% Deviation mass fluxes
F1_pert = work.F1 - m1bar.*w_bar_star;
F2_pert = work.F2 - m2bar.*w_bar_star;

% Mean ascent contribution
QF_mean = rhow.*w_bar_star.*q_bar_star;

% Resolved contributions in fluids 1 and 2
QF_res1 = F1_pert.*(q1 - q_bar_star);
QF_res2 = F2_pert.*(q2 - q_bar_star);

% Subfilter-scale contributions
QF_sf1 = Dq1;
QF_sf2 = Dq2;
QF_sf1_bc = work.Dq1bc;
QF_sf2_bc = work.Dq2bc;
QF_sf1_ed = work.Dq1ed;
QF_sf2_ed = work.Dq2ed;

% Total
QF_tot = QF_mean + QF_res1 + QF_res2 + weight_to_w(grid,QF_sf1 + QF_sf2);

% And plot

subplot(1,3,1)
plot(QF_tot,zunitsw,'k',...
     QF_mean,zunitsw,'r')
ylim([0,zplottop])
title('Tot and expansion')
xlabel('w q')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(1,3,2)
plot(QF_res1,zunitsw,'b',...
     QF_res2,zunitsw,'r')
ylim([0,zplottop])
title('Resolved')
xlabel('w q')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(1,3,3)
plot(QF_sf1_ed,zunitsp,'b',...
     QF_sf2_ed,zunitsp,'r',...
     QF_sf1_bc,zunitsp,'b--',...
     QF_sf2_bc,zunitsp,'r--')
ylim([0,zplottop])
title('Subfilter')
xlabel('w q')
ylabel(labelz)
set(gca,'FontSize',fs)

