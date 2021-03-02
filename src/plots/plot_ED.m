% Plot actual and candidate entrainment and detrainment profiles


% ----------

% Figure 18: actual and candidate entrainment and detrainment profiles

figure(18)

% subplot(1,4,1)
% plot(relabel.M21,zunitsp,'r',...
%     -relabel.M12,zunitsp,'b',...
%      relabel.trialM21_buoy,zunitsp,'r--',...
%     -relabel.trialM12_buoy,zunitsp,'b--')
% ylim([0,zplottop])
% title('E and D')
% xlabel('dm2/dt')
% ylabel(labelz)
% legend('E','D','Eb','Db','Location','NorthEast')
% set(gca,'FontSize',fs)

subplot(2,4,1)
plot(relabel.M21,zunitsp,'r',...
    -relabel.M12,zunitsp,'b',...
    -relabel.M12_mix ,zunitsp,'k--',...
    -relabel.M12_sort,zunitsp,'k')
ylim([0,zplottop])
title('E and D')
xlabel('dm2/dt')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(2,4,5)
plot(relabel.f_sort,zunitsw,'k')
ylim([0,zplottop])
title('f_{sort}')
xlabel('f_{sort}')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(2,4,2)
plot(relabel.what21,zunitsw,'k')
ylim([0,zplottop])
title('$$\hat{w}_{21}$$','interpreter','latex')
xlabel('what21')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(2,4,6)
plot(relabel.what12_mix,  zunitsw,'k--',...
     relabel.what12_sort ,zunitsw,'k',...
     relabel.what12_blend,zunitsw,'b')
ylim([0,zplottop])
title('$$\hat{w}_{12}$$','interpreter','latex')
xlabel('what12')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(2,4,3)
plot(relabel.etahat21,zunitsw,'k')
ylim([0,zplottop])
title('$$\hat{\eta}_{21}$$','interpreter','latex')
xlabel('etahat21')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(2,4,7)
plot(relabel.etahat12_mix,  zunitsw,'k--',...
     relabel.etahat12_sort ,zunitsw,'k',...
     relabel.etahat12_blend,zunitsw,'b')
ylim([0,zplottop])
title('$$\hat{\eta}_{12}$$','interpreter','latex')
xlabel('etahat12')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(2,4,4)
plot(relabel.qhat21,zunitsw,'k')
ylim([0,zplottop])
title('$$\hat{q}_{21}$$','interpreter','latex')
xlabel('qhat21')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(2,4,8)
plot(relabel.qhat12_mix,  zunitsw,'k--',...
     relabel.qhat12_sort ,zunitsw,'k',...
     relabel.qhat12_blend,zunitsw,'b')
ylim([0,zplottop])
title('$$\hat{q}_{12}$$','interpreter','latex')
xlabel('qhat12')
ylabel(labelz)
set(gca,'FontSize',fs)




% pause
