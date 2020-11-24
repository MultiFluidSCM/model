% Plot diffusive and buoyancy correlation eta flux


% ----------

% Figure 25: diffusive and buoyancy corelation eta flux

figure(25)

subplot(1,1,1)
plot(work.Deta1,zunitsp,'b',...
     work.Deta2,zunitsp,'r',...
     work.Deta1 + work.Deta2,zunitsp,'k',...
     work.Deta1bc,zunitsp,'b--',...
     work.Deta2bc,zunitsp,'r--',...
     work.Deta1bc + work.Deta2bc,zunitsp,'k--')
ylim([0,zplottop])
title('SG eta flux')
xlabel('w \eta')
ylabel(labelz)
set(gca,'FontSize',fs)

