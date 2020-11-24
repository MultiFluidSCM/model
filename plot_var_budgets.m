% Plot variance budgets


% ----------

% Figure 26: Budgets for eta and q variances

figure(26)

subplot(2,2,1)
plot(tend.fluid(1).mvareta.diffuse,zunitsw,'g',...
     tend.fluid(1).mvareta.diffent,zunitsw,'g--',...
     tend.fluid(1).mvareta.dissn  ,zunitsw,'r',...
     tend.fluid(1).mvareta.relabel,zunitsw,'b',...
     tend.fluid(1).mvareta.tot    ,zunitsw,'ko')
ylim([0,zplottop])
title('var eta1 budget')
xlabel('dmVe/dt')
ylabel(labelz)
legend('Diff','Dife','Diss','R','Tot','Location','NorthEast')
set(gca,'FontSize',fs)

subplot(2,2,2)
plot(tend.fluid(2).mvareta.diffuse,zunitsw,'g',...
     tend.fluid(2).mvareta.diffent,zunitsw,'g--',...
     tend.fluid(2).mvareta.dissn  ,zunitsw,'r',...
     tend.fluid(2).mvareta.relabel,zunitsw,'b',...
     tend.fluid(2).mvareta.tot    ,zunitsw,'ko')
ylim([0,zplottop])
title('var eta2 budget')
xlabel('dmVe/dt')
ylabel(labelz)
legend('Diff','Dife','Diss','R','Tot','Location','NorthEast')
set(gca,'FontSize',fs)

subplot(2,2,3)
plot(tend.fluid(1).mvarq.diffuse,zunitsw,'g',...
     tend.fluid(1).mvarq.diffent,zunitsw,'g--',...
     tend.fluid(1).mvarq.dissn  ,zunitsw,'r',...
     tend.fluid(1).mvarq.relabel,zunitsw,'b',...
     tend.fluid(1).mvarq.tot    ,zunitsw,'ko')
ylim([0,zplottop])
title('var q1 budget')
xlabel('dmVq/dt')
ylabel(labelz)
legend('Diff','Dife','Diss','R','Tot','Location','NorthEast')
set(gca,'FontSize',fs)

subplot(2,2,4)
plot(tend.fluid(2).mvarq.diffuse,zunitsw,'g',...
     tend.fluid(2).mvarq.diffent,zunitsw,'g--',...
     tend.fluid(2).mvarq.dissn  ,zunitsw,'r',...
     tend.fluid(2).mvarq.relabel,zunitsw,'b',...
     tend.fluid(2).mvarq.tot    ,zunitsw,'ko')
ylim([0,zplottop])
title('var q2 budget')
xlabel('dmVq/dt')
ylabel(labelz)
legend('Diff','Dife','Diss','R','Tot','Location','NorthEast')
set(gca,'FontSize',fs)



