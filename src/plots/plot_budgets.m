% Plot updraft budgets


% ----------

% Figure 3: Updraft budgets for mass, w, eta, q
% Also TKE budgets

% Updraft mass budget
if settings.switches.plot_budgets_mass
    if plottype == 0
        figure(3)
        subplot(2,4,1)
        plot(tend.fluid(2).m.transport,zunitsp,'g',...
             relabel.M21,zunitsp,'r',...
            -relabel.M12,zunitsp,'b',...
             tend.fluid(2).m.tot,zunitsp,'k')
        ylim([0,zplottop])
        title('m_2 budget')
        xlabel('dm2/dt')
        ylabel(labelz)
        legend('Tr','E','D','To','Location','NorthEast')
        set(gca,'FontSize',fs)
    else
        figure(5)
        subplot(1,1,1)
        plot(tend.fluid(2).m.transport,zunitsp,'g',...
             relabel.M21,zunitsp,'r',...
            -relabel.M12,zunitsp,'b',...
             tend.fluid(2).m.tot,zunitsp,'k')
        ylim([0,zplottop])
        title('m_2 budget')
        xlabel('dm2/dt')
        ylabel(labelz)
        legend('Tr','E','D','To','Location','NorthEast')
        set(gca,'FontSize',fs)
        pause
    end
end

% Needed for several diagnostics below
m1transbar = weight_to_w(grid,tend.fluid(1).m.transport);
m1totbar   = weight_to_w(grid,tend.fluid(1).m.tot);
m2transbar = weight_to_w(grid,tend.fluid(2).m.transport);
m2totbar   = weight_to_w(grid,tend.fluid(2).m.tot);

% Updraft vertical velocity budget (advective form)
if settings.switches.plot_budgets_w
    if plottype == 0
        figure(3)
        subplot(2,4,5)
        plot(budgets.w2.transport,zunitsw,'g',...
             budgets.w2.diffuse  ,zunitsw,'g--',...
             budgets.w2.entrain  ,zunitsw,'r',...
             budgets.w2.detrain  ,zunitsw,'b',...
             budgets.w2.pgterm   ,zunitsw,'b--',...
             budgets.w2.drag     ,zunitsw,'r--',...
             budgets.w2.wfix     ,zunitsw,'k--',...
             budgets.w2.tot      ,zunitsw,'k')
        ylim([0,zplottop])
        title('w_2 budget')
        xlabel('dw_2/dt')
        ylabel(labelz)
        legend('Tr','K','E','D','PG','P','X','To','Location','NorthEast')
        set(gca,'FontSize',fs)
    else
        figure(5)
        subplot(1,1,1)
        plot(budgets.w2.transport,zunitsw,'g',...
             budgets.w2.diffuse  ,zunitsw,'g--',...
             budgets.w2.entrain  ,zunitsw,'r',...
             budgets.w2.detrain  ,zunitsw,'b',...
             budgets.w2.pgterm   ,zunitsw,'b--',...
             budgets.w2.drag     ,zunitsw,'r--',...
             budgets.w2.wfix     ,zunitsw,'k--',...
             budgets.w2.tot      ,zunitsw,'k')
        ylim([0,zplottop])
        title('w_2 budget')
        xlabel('dw_2/dt')
        ylabel(labelz)
        legend('Tr','K','E','D','PG','P','X','To','Location','NorthEast')
        set(gca,'FontSize',fs)
        pause
    end
end

if settings.switches.plot_budgets_entropy
    % Updraft entropy budget (advective form)
    if plottype == 0
        figure(3)
        subplot(2,4,6)
        plot(budgets.eta2.transport,zunitsw,'g',...
             budgets.eta2.diffuse  ,zunitsw,'g--',...
             budgets.eta2.entrain  ,zunitsw,'r',...
             budgets.eta2.detrain  ,zunitsw,'b',...
             budgets.eta2.bflux    ,zunitsw,'b--',...
             budgets.eta2.dissn    ,zunitsw,'r--',...
             budgets.eta2.wfix     ,zunitsw,'k--',...
             budgets.eta2.tot      ,zunitsw,'k')
        ylim([0,zplottop])
        title('\eta_2 budget')
        xlabel('deta_2/dt')
        ylabel(labelz)
        legend('Tr','K','E','D','B','Q','X','To','Location','NorthEast')
        set(gca,'FontSize',14)
    else
        figure(5)
        subplot(1,1,1)
        plot(budgets.eta2.transport,zunitsw,'g',...
             budgets.eta2.diffuse  ,zunitsw,'g--',...
             budgets.eta2.entrain  ,zunitsw,'r',...
             budgets.eta2.detrain  ,zunitsw,'b',...
             budgets.eta2.bflux    ,zunitsw,'b--',...
             budgets.eta2.dissn    ,zunitsw,'r--',...
             budgets.eta2.wfix     ,zunitsw,'k--',...
             budgets.eta2.tot      ,zunitsw,'k')
        ylim([0,zplottop])
        title('\eta_2 budget')
        xlabel('deta_2/dt')
        ylabel(labelz)
        legend('Tr','K','E','D','B','Q','X','To','Location','NorthEast')
        set(gca,'FontSize',fs)
        pause
    end


    % Downdraft entropy budget (advective form)
    if plottype == 0
        figure(3)
        subplot(2,4,2)
        plot(budgets.eta1.transport,zunitsw,'g',...
             budgets.eta1.diffuse  ,zunitsw,'g--',...
             budgets.eta1.entrain  ,zunitsw,'r',...
             budgets.eta1.detrain  ,zunitsw,'b',...
             budgets.eta1.bflux    ,zunitsw,'b--',...
             budgets.eta1.dissn    ,zunitsw,'r--',...
             budgets.eta1.wfix     ,zunitsw,'k--',...
             budgets.eta1.tot      ,zunitsw,'k')
        ylim([0,zplottop])
        title('\eta_1 budget')
        xlabel('deta_1/dt')
        ylabel(labelz)
        legend('Tr','K','E','D','B','Q','X','To','Location','NorthEast')
        set(gca,'FontSize',fs)
    else
        figure(5)
        subplot(1,1,1)
        plot(budgets.eta1.transport,zunitsw,'g',...
             budgets.eta1.diffuse  ,zunitsw,'g--',...
             budgets.eta1.entrain  ,zunitsw,'r',...
             budgets.eta1.detrain  ,zunitsw,'b',...
             budgets.eta1.bflux    ,zunitsw,'b--',...
             budgets.eta1.dissn    ,zunitsw,'r--',...
             budgets.eta1.wfix     ,zunitsw,'k--',...
             budgets.eta1.tot      ,zunitsw,'k')
        ylim([0,zplottop])
        title('\eta_1 budget')
        xlabel('deta_1/dt')
        ylabel(labelz)
        legend('Tr','K','E','D','B','Q','X','To','Location','NorthEast')
        set(gca,'FontSize',fs)
        pause
    end
end

if settings.switches.plot_budgets_water
    % Updraft water budget (advective form)
    if plottype == 0
        figure(3)
        subplot(2,4,7)
        plot(budgets.q2.transport,zunitsw,'g',...
             budgets.q2.diffuse  ,zunitsw,'g--',...
             budgets.q2.entrain  ,zunitsw,'r',...
             budgets.q2.detrain  ,zunitsw,'b',...
             budgets.q2.wfix     ,zunitsw,'k--',...
             budgets.q2.tot      ,zunitsw,'k')
        ylim([0,zplottop])
        title('q_2 budget')
        xlabel('dq_2/dt')
        ylabel(labelz)
        legend('Tr','K','E','D','X','To','Location','NorthEast')
        set(gca,'FontSize',fs)
    else
        figure(5)
        subplot(1,1,1)
        plot(budgets.q2.transport,zunitsw,'g',...
             budgets.q2.diffuse  ,zunitsw,'g--',...
             budgets.q2.entrain  ,zunitsw,'r',...
             budgets.q2.detrain  ,zunitsw,'b',...
             budgets.q2.wfix     ,zunitsw,'k--',...
             budgets.q2.tot      ,zunitsw,'k')
        ylim([0,zplottop])
        title('q_2 budget')
        xlabel('dq_2/dt')
        ylabel(labelz)
        legend('Tr','K','E','D','X','To','Location','NorthEast')
        set(gca,'FontSize',fs)
        pause
    end


    % Downdraft water budget (advective form)
    if plottype == 0
        figure(3)
        subplot(2,4,3)
        plot(budgets.q1.transport,zunitsw,'g',...
             budgets.q1.diffuse  ,zunitsw,'g--',...
             budgets.q1.entrain  ,zunitsw,'r',...
             budgets.q1.detrain  ,zunitsw,'b',...
             budgets.q1.wfix     ,zunitsw,'k--',...
             budgets.q1.tot      ,zunitsw,'k')
        ylim([0,zplottop])
        title('q_1 budget')
        xlabel('dq_1/dt')
        ylabel(labelz)
        legend('Tr','K','E','D','X','To','Location','NorthEast')
        set(gca,'FontSize',fs)
    else
        figure(5)
        subplot(1,1,1)
        plot(budgets.q1.transport,zunitsw,'g',...
             budgets.q1.diffuse  ,zunitsw,'g--',...
             budgets.q1.entrain  ,zunitsw,'r',...
             budgets.q1.detrain  ,zunitsw,'b',...
             budgets.q1.wfix     ,zunitsw,'k--',...
             budgets.q1.tot      ,zunitsw,'k')
        ylim([0,zplottop])
        title('q_1 budget')
        xlabel('dq_1/dt')
        ylabel(labelz)
        legend('Tr','K','E','D','X','To','Location','NorthEast')
        set(gca,'FontSize',fs)
        pause
    end
end

% TKE budgets (flux form)
if settings.switches.plot_budgets_tke
    if plottype == 0
        figure(3)
        subplot(2,4,4)
        plot(tend.fluid(1).mtke.transport,zunitsp,'g',...
             tend.fluid(1).mtke.relabel,zunitsp,'r',...
             tend.fluid(1).mtke.shear,zunitsp,'g--',...
             tend.fluid(1).mtke.bflux,zunitsp,'b',...
             tend.fluid(1).mtke.diffuse + tend.fluid(1).mtke.diffent,zunitsp,'b:',...
             tend.fluid(1).mtke.drag,zunitsp,'k--',...
             tend.fluid(1).mtke.dissn,zunitsp,'r--',...
             tend.fluid(1).mtke.tot,zunitsp,'k')
        ylim([0,zplottop])
        title('TKE_1 budget')
        xlabel('d/dt')
        ylabel(labelz)
        legend('Tr','E-D','S','B','K','P','D','To','Location','NorthEast')
        set(gca,'FontSize',fs)
    else
        figure(5)
        subplot(1,1,1)
        plot(tend.fluid(1).mtke.transport,zunitsp,'g',...
             tend.fluid(1).mtke.relabel,zunitsp,'r',...
             tend.fluid(1).mtke.shear,zunitsp,'g--',...
             tend.fluid(1).mtke.bflux,zunitsp,'b',...
             tend.fluid(1).mtke.diffuse + tend.fluid(1).mtke.diffent,zunitsp,'b:',...
             tend.fluid(1).mtke.drag,zunitsp,'k--',...
             tend.fluid(1).mtke.dissn,zunitsp,'r--',...
             tend.fluid(1).mtke.tot,zunitsp,'k')
        ylim([0,zplottop])
        title('TKE_1 budget')
        xlabel('d/dt')
        ylabel(labelz)
        legend('Tr','E-D','S','B','K','P','D','To','Location','NorthEast')
        set(gca,'FontSize',fs)
        pause
    end
    if plottype == 0
        figure(3)
        subplot(2,4,8)
        plot(tend.fluid(2).mtke.transport,zunitsp,'g',...
             tend.fluid(2).mtke.relabel,zunitsp,'r',...
             tend.fluid(2).mtke.shear,zunitsp,'g--',...
             tend.fluid(2).mtke.bflux,zunitsp,'b',...
             tend.fluid(2).mtke.diffuse + tend.fluid(2).mtke.diffent,zunitsp,'b:',...
             tend.fluid(2).mtke.drag,zunitsp,'k--',...
             tend.fluid(2).mtke.dissn,zunitsp,'r--',...
             tend.fluid(2).mtke.tot,zunitsp,'k')
        ylim([0,zplottop])
        title('TKE_2 budget')
        xlabel('d/dt')
        ylabel(labelz)
        legend('Tr','E-D','S','B','K','P','D','To','Location','NorthEast')
        set(gca,'FontSize',fs)
    else
        figure(5)
        subplot(1,1,1)
        plot(tend.fluid(2).mtke.transport,zunitsp,'g',...
             tend.fluid(2).mtke.relabel,zunitsp,'r',...
             tend.fluid(2).mtke.shear,zunitsp,'g--',...
             tend.fluid(2).mtke.bflux,zunitsp,'b',...
             tend.fluid(2).mtke.diffuse + tend.fluid(2).mtke.diffent,zunitsp,'b:',...
             tend.fluid(2).mtke.drag,zunitsp,'k--',...
             tend.fluid(2).mtke.dissn,zunitsp,'r--',...
             tend.fluid(2).mtke.tot,zunitsp,'k')
        ylim([0,zplottop])
        title('TKE_2 budget')
        xlabel('d/dt')
        ylabel(labelz)
        legend('Tr','E-D','S','B','K','P','D','To','Location','NorthEast')
        set(gca,'FontSize',fs)
        pause
    end
end

% Entrainment budget
if settings.switches.plot_budgets_transfers
    fig = figure(5);
    set(gcf,'Position',[76 166 1257 624])
    
    relabelMax = max(max(relabel.M21, relabel.M12));
    
    subplot(1,3,1)
    plot(relabel.M21_instab,zunitsp,'r',...
         relabel.M21_sort,  zunitsp,'b',...
         relabel.M21_dwdz,  zunitsp,'c',...
         relabel.M21_mix,   zunitsp,'m',...
         relabel.M21,       zunitsp,'k--')
    xlim([0,relabelMax])
    ylim([0,zplottop])
    title('Entrainment')
    xlabel('dm2/dt')
    ylabel(labelz)
    legend('Instab','Sort','dw/dz','Mix','Total','Location','NorthEast')
    set(gca,'FontSize',fs)
    
    subplot(1,3,2)
    plot(relabel.M12_instab,zunitsp,'r',...
         relabel.M12_sort,  zunitsp,'b',...
         relabel.M12_dwdz,  zunitsp,'c',...
         relabel.M12_mix,   zunitsp,'m',...
         relabel.M12,       zunitsp,'k--')
    xlim([0,relabelMax])
    ylim([0,zplottop])
    title('Detrainment')
    xlabel('dm2/dt')
    ylabel(labelz)
    legend('Instab','Sort','dw/dz','Mix','Total','Location','NorthEast')
    set(gca,'FontSize',fs)
    
    if isfield(ts, 'zstar')
        subplot(1,3,3)
        % Cloud top and base
        % Tidy cases with cloud base at lid
        SCM_zstar = ts.zstar*1e-3;
        SCM_zcbase = ts.zcbaseSG*1e-3;
        SCM_zctop = ts.zctopSG*1e-3;
        SCM_time_ser_hours = ts.time_high_res/3600.;
        plot(SCM_time_ser_hours,SCM_zstar ,'b--',...
             SCM_time_ser_hours,SCM_zctop ,'r--',...
             SCM_time_ser_hours,SCM_zcbase,'k--')
        xlabel('Time','fontsize',fs)
        ylabel('Cld base/top (km)','fontsize',fs)
        xlim([0,15]);
        ylim([0,zplottop]);
        %title('Cloud base/top','fontsize',fs)
        sgtitle([num2str(ts.time(end)/3600,'%.2f'),' hours'],'fontsize',1.5*fs)
        set(gca,'fontsize',fs,'XTick',[1:3:14])
    end
    
    saveas(...
        fig,...
        fullfile(...
            settings.folders.images,...
            join(["profiles_budget_transfers_",num2str(plotstep),".png"], "")...
        )...
    );
end
