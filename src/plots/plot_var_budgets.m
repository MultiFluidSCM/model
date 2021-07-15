% Plot variance budgets


% ----------

% Figure 26: Budgets for eta and q variances and covariances

figure(26)

subplot(3,2,1)
plot(tend.fluid(1).mvareta.diffuse,zunitsp,'g',...
     tend.fluid(1).mvareta.diffent,zunitsp,'g--',...
     tend.fluid(1).mvareta.dissn  ,zunitsp,'r',...
     tend.fluid(1).mvareta.relabel,zunitsp,'b',...
     tend.fluid(1).mvareta.tot    ,zunitsp,'ko')
ylim([0,zplottop])
title('var eta1 budget')
xlabel('dmVe/dt')
ylabel(labelz)
legend('Diff','Dife','Diss','R','Tot','Location','NorthEast')
set(gca,'FontSize',fs)

subplot(3,2,2)
plot(tend.fluid(2).mvareta.diffuse,zunitsp,'g',...
     tend.fluid(2).mvareta.diffent,zunitsp,'g--',...
     tend.fluid(2).mvareta.dissn  ,zunitsp,'r',...
     tend.fluid(2).mvareta.relabel,zunitsp,'b',...
     tend.fluid(2).mvareta.tot    ,zunitsp,'ko')
ylim([0,zplottop])
title('var eta2 budget')
xlabel('dmVe/dt')
ylabel(labelz)
legend('Diff','Dife','Diss','R','Tot','Location','NorthEast')
set(gca,'FontSize',fs)

subplot(3,2,3)
plot(tend.fluid(1).mvarq.diffuse,zunitsp,'g',...
     tend.fluid(1).mvarq.diffent,zunitsp,'g--',...
     tend.fluid(1).mvarq.dissn  ,zunitsp,'r',...
     tend.fluid(1).mvarq.relabel,zunitsp,'b',...
     tend.fluid(1).mvarq.tot    ,zunitsp,'ko')
ylim([0,zplottop])
title('var q1 budget')
xlabel('dmVq/dt')
ylabel(labelz)
legend('Diff','Dife','Diss','R','Tot','Location','NorthEast')
set(gca,'FontSize',fs)

subplot(3,2,4)
plot(tend.fluid(2).mvarq.diffuse,zunitsp,'g',...
     tend.fluid(2).mvarq.diffent,zunitsp,'g--',...
     tend.fluid(2).mvarq.dissn  ,zunitsp,'r',...
     tend.fluid(2).mvarq.relabel,zunitsp,'b',...
     tend.fluid(2).mvarq.tot    ,zunitsp,'ko')
ylim([0,zplottop])
title('var q2 budget')
xlabel('dmVq/dt')
ylabel(labelz)
legend('Diff','Dife','Diss','R','Tot','Location','NorthEast')
set(gca,'FontSize',fs)

subplot(3,2,5)
plot(tend.fluid(1).mcovaretaq.diffuse,zunitsp,'g',...
     tend.fluid(1).mcovaretaq.diffent,zunitsp,'g--',...
     tend.fluid(1).mcovaretaq.dissn  ,zunitsp,'r',...
     tend.fluid(1).mcovaretaq.relabel,zunitsp,'b',...
     tend.fluid(1).mcovaretaq.tot    ,zunitsp,'ko')
ylim([0,zplottop])
title('covar 1 budget')
xlabel('dmC/dt')
ylabel(labelz)
legend('Diff','Dife','Diss','R','Tot','Location','NorthEast')
set(gca,'FontSize',fs)

subplot(3,2,6)
plot(tend.fluid(2).mcovaretaq.diffuse,zunitsp,'g',...
     tend.fluid(2).mcovaretaq.diffent,zunitsp,'g--',...
     tend.fluid(2).mcovaretaq.dissn  ,zunitsp,'r',...
     tend.fluid(2).mcovaretaq.relabel,zunitsp,'b',...
     tend.fluid(2).mcovaretaq.tot    ,zunitsp,'ko')
ylim([0,zplottop])
title('covar 2 budget')
xlabel('dmC/dt')
ylabel(labelz)
legend('Diff','Dife','Diss','R','Tot','Location','NorthEast')
set(gca,'FontSize',fs)

