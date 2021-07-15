% Plot covariances

% ----------

% Figure 34: covariances

% The covariances should be bounded by +/- sqrt(vareta*varq),
% so plot these bounds for reference

figure(34)

covarmax1 = sqrt(state_new.fluid(1).covaretaq.*state_new.fluid(1).covaretaq);
covarmax2 = sqrt(state_new.fluid(2).covaretaq.*state_new.fluid(2).covaretaq);

% And plot

subplot(1,2,1)
plot(-covarmax1,                   zunitsp,'r--',...
      covarmax1,                   zunitsp,'r--',...
      state_new.fluid(1).covaretaq,zunitsp,'k')
ylim([0,zplottop])
title('covar 1')
xlabel('\eta q')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(1,2,2)
plot(-covarmax2,                   zunitsp,'r--',...
      covarmax2,                   zunitsp,'r--',...
      state_new.fluid(2).covaretaq,zunitsp,'k')
ylim([0,zplottop])
title('covar 2')
xlabel('\eta q')
ylabel(labelz)
set(gca,'FontSize',fs)
