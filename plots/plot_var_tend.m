% Plot eta and q variance budgets


% ----------

% Figure 23: Budgets for eta and q variances

% Updraft eta variance budget
if plottype == 0
    figure(23)
    subplot(2,2,1)
    plot(tend.fluid(2).mvareta.diffuse,zunitsp,'g',...
         tend.fluid(2).mvareta.diffent,zunitsp,'g--',...
         tend.fluid(2).mvareta.dissn,zunitsp,'r',...
         tend.fluid(2).mvareta.relabel,zunitsp,'b',...
         tend.fluid(2).mvareta.tot,zunitsp,'k')
    ylim([0,zplottop])
    title('var-eta_2 budget')
    xlabel('m*d/dt')
    ylabel(labelz)
    legend('K','Ke','D','R','Tot','Location','NorthEast')
    set(gca,'FontSize',fs)
else
    figure(5)
    subplot(1,1,1)
    plot(tend.fluid(2).mvareta.diffuse,zunitsp,'g',...
         tend.fluid(2).mvareta.diffent,zunitsp,'g--',...
         tend.fluid(2).mvareta.dissn,zunitsp,'r',...
         tend.fluid(2).mvareta.relabel,zunitsp,'b',...
         tend.fluid(2).mvareta.tot,zunitsp,'k')
    ylim([0,zplottop])
    title('var-eta_2 budget')
    xlabel('m*d/dt')
    ylabel(labelz)
    legend('K','Ke','D','R','Tot','Location','NorthEast')
    set(gca,'FontSize',fs)
    pause
end


% Downdraft eta variance budget
if plottype == 0
    figure(23)
    subplot(2,2,3)
    plot(tend.fluid(1).mvareta.diffuse,zunitsp,'g',...
         tend.fluid(1).mvareta.diffent,zunitsp,'g--',...
         tend.fluid(1).mvareta.dissn,zunitsp,'r',...
         tend.fluid(1).mvareta.relabel,zunitsp,'b',...
         tend.fluid(1).mvareta.tot,zunitsp,'k')
    ylim([0,zplottop])
    title('var-eta_1 budget')
    xlabel('m*d/dt')
    ylabel(labelz)
    legend('K','Ke','D','R','Tot','Location','NorthEast')
    set(gca,'FontSize',fs)
else
    figure(5)
    subplot(1,1,1)
    plot(tend.fluid(1).mvareta.diffuse,zunitsp,'g',...
         tend.fluid(1).mvareta.diffent,zunitsp,'g--',...
         tend.fluid(1).mvareta.dissn,zunitsp,'r',...
         tend.fluid(1).mvareta.relabel,zunitsp,'b',...
         tend.fluid(1).mvareta.tot,zunitsp,'k')
    ylim([0,zplottop])
    title('var-eta_1 budget')
    xlabel('m*d/dt')
    ylabel(labelz)
    legend('K','Ke','D','R','Tot','Location','NorthEast')
    set(gca,'FontSize',fs)
    pause
end

% Updraft q variance budget
if plottype == 0
    figure(23)
    subplot(2,2,2)
    plot(tend.fluid(2).mvarq.diffuse,zunitsp,'g',...
         tend.fluid(2).mvarq.diffent,zunitsp,'g--',...
         tend.fluid(2).mvarq.dissn,zunitsp,'r',...
         tend.fluid(2).mvarq.relabel,zunitsp,'b',...
         tend.fluid(2).mvarq.tot,zunitsp,'k')
    ylim([0,zplottop])
    title('var-q_2 budget')
    xlabel('m*d/dt')
    ylabel(labelz)
    legend('K','Ke','D','R','Tot','Location','NorthEast')
    set(gca,'FontSize',fs)
else
    figure(5)
    subplot(1,1,1)
    plot(tend.fluid(2).mvarq.diffuse,zunitsp,'g',...
         tend.fluid(2).mvarq.diffent,zunitsp,'g--',...
         tend.fluid(2).mvarq.dissn,zunitsp,'r',...
         tend.fluid(2).mvarq.relabel,zunitsp,'b',...
         tend.fluid(2).mvarq.tot,zunitsp,'k')
    ylim([0,zplottop])
    title('var-q_2 budget')
    xlabel('m*d/dt')
    ylabel(labelz)
    legend('K','Ke','D','R','Tot','Location','NorthEast')
    set(gca,'FontSize',fs)
    pause
end

% Downdraft q variance budget
if plottype == 0
    figure(23)
    subplot(2,2,4)
    plot(tend.fluid(1).mvarq.diffuse,zunitsp,'g',...
         tend.fluid(1).mvarq.diffent,zunitsp,'g--',...
         tend.fluid(1).mvarq.dissn,zunitsp,'r',...
         tend.fluid(1).mvarq.relabel,zunitsp,'b',...
         tend.fluid(1).mvarq.tot,zunitsp,'k')
    ylim([0,zplottop])
    title('var-q_1 budget')
    xlabel('m*d/dt')
    ylabel(labelz)
    legend('K','Ke','D','R','Tot','Location','NorthEast')
    set(gca,'FontSize',fs)
else
    figure(5)
    subplot(1,1,1)
    plot(tend.fluid(1).mvarq.diffuse,zunitsp,'g',...
         tend.fluid(1).mvarq.diffent,zunitsp,'g--',...
         tend.fluid(1).mvarq.dissn,zunitsp,'r',...
         tend.fluid(1).mvarq.relabel,zunitsp,'b',...
         tend.fluid(1).mvarq.tot,zunitsp,'k')
    ylim([0,zplottop])
    title('var-q_1 budget')
    xlabel('m*d/dt')
    ylabel(labelz)
    legend('K','Ke','D','R','Tot','Location','NorthEast')
    set(gca,'FontSize',fs)
    pause
end

