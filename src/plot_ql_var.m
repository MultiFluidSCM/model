% Plot subgrid variances related to liquid water

% ----------

% Figure 28: Variances related to liquid water

figure(28)

subplot(2,3,1)
plot(state_new.fluid(1).varq,zunitsw,'b',...
     state_new.fluid(2).varq,zunitsw,'r')
ylim([0,zplottop])
title('Var q')
xlabel('Var q')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(2,3,4)
plot(state_new.fluid(1).vareta,zunitsw,'b',...
     state_new.fluid(2).vareta,zunitsw,'r')
ylim([0,zplottop])
title('Var eta')
xlabel('Var eta')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(2,3,2)
plot(eos.Covarqlq1,zunitsw,'b',...
     eos.Covarqlq2,zunitsw,'r')
ylim([0,zplottop])
title('Covar ql-q')
xlabel('Covar ql-q')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(2,3,5)
plot(eos.Covarqleta1,zunitsw,'b',...
     eos.Covarqleta2,zunitsw,'r')
ylim([0,zplottop])
title('Covar ql-eta')
xlabel('Covar ql-eta')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(2,3,3)
plot(eos.Varql1,zunitsw,'b',...
     eos.Varql2,zunitsw,'r')
ylim([0,zplottop])
title('Var ql')
xlabel('Var ql')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(2,3,6)
plot(eos.Varrhow1,zunitsw,'b',...
     eos.Varrhow2,zunitsw,'r',...
     eos.ccc1,zunitsw,'b--',...
     eos.ccc2,zunitsw,'k--',...
     eos.ccc3,zunitsw,'r--')
ylim([0,zplottop])
title('Var rho')
xlabel('Var rho')
ylabel(labelz)
set(gca,'FontSize',fs)


% pause
