% Plot a variety of turbulence related diagnostics


% ----------

% Figure 4: TKE, turbulence length and time scales
% and implied diffusion coefficient


% TKE
if plottype == 0
    figure(4)
    subplot(2,4,1)
    plot(tke1,zunitsp,'k')
    % plot(tke1,zunitsp,'k',[scales.tke_thresh,scales.tke_thresh],[0,zplottop],'k--')
    ylim([0,zplottop])
    title('TKE_1')
    xlabel('TKE_1')
    ylabel(labelz)
    set(gca,'FontSize',fs)
else
    figure(5)
    subplot(1,1,1)
    plot(tke1,zunitsp,'k')
    ylim([0,zplottop])
    title('TKE_1')
    xlabel('TKE_1')
    ylabel(labelz)
    set(gca,'FontSize',fs)
    pause
end
if plottype == 0
    figure(4)
    subplot(2,4,5)
    plot(tke2,zunitsp,'k')
    ylim([0,zplottop])
    title('TKE_2')
    xlabel('TKE_2')
    ylabel(labelz)
    set(gca,'FontSize',fs)
else
    figure(5)
    subplot(1,1,1)
    plot(tke2,zunitsp,'k')
    ylim([0,zplottop])
    title('TKE_2')
    xlabel('TKE_2')
    ylabel(labelz)
    set(gca,'FontSize',fs)
    pause
end


% Turbulence and plume length scales
if plottype == 0
    figure(4)
    subplot(2,4,2)
    plot(scales.L_turb1,zunitsp,'k')
    ylim([0,zplottop])
    title('L_1')
    xlabel('L_1')
    ylabel(labelz)
    set(gca,'FontSize',fs)
else
    figure(5)
    subplot(1,1,1)
    plot(scales.L_turb1,zunitsp,'k')
    ylim([0,zplottop])
    title('L_1')
    xlabel('L_1')
    ylabel(labelz)
    set(gca,'FontSize',fs)
    pause
end
if plottype == 0
    figure(4)
    subplot(2,4,6)
    plot(scales.L_turb2,zunitsp,'k',scales.L_plume,zunitsp,'k--')
    ylim([0,zplottop])
    title('L_2, L_{plume}')
    xlabel('L_2')
    ylabel(labelz)
    set(gca,'FontSize',fs)
else
    figure(5)
    subplot(1,1,1)
    plot(scales.L_turb2,zunitsp,'k',scales.L_plume,zunitsp,'k--')
    ylim([0,zplottop])
    title('L_2, L_{plume}')
    xlabel('L_2')
    ylabel(labelz)
    set(gca,'FontSize',fs)
    pause
end

% dw2dz = (w2(2:nzp) - w2(1:nz))./grid.dzp;
% Turbulence time scale
if plottype == 0
    figure(4)
    subplot(2,4,3)
    plot(scales.T_turb1,zunitsp,'k')
    ylim([0,zplottop])
    title('T_1')
    xlabel('T_1')
    ylabel(labelz)
    set(gca,'FontSize',fs)
else
    figure(5)
    subplot(1,1,1)
    plot(scales.T_turb1,zunitsp,'k')
    ylim([0,zplottop])
    title('T_1')
    xlabel('T_1')
    ylabel(labelz)
    set(gca,'FontSize',fs)
    pause
end
if plottype == 0
    figure(4)
    subplot(2,4,7)
    plot(scales.T_turb2,zunitsp,'k')
    %plot(1./scales.T_turb2,zunitsp,'k',-dw2dz,zunitsp,'r',relabel.rate_mix,zunitsp,'r--')
    ylim([0,zplottop])
    title('T_2')
    xlabel('T_2')
    ylabel(labelz)
    set(gca,'FontSize',fs)
else
    figure(5)
    subplot(1,1,1)
    plot(scales.T_turb2,zunitsp,'k')
    ylim([0,zplottop])
    title('T_2')
    xlabel('T_2')
    ylabel(labelz)
    set(gca,'FontSize',fs)
    pause
end


% Implied diffusivity cf actual diffusivity
if plottype == 0
    figure(4)
    subplot(2,4,4)
    plot(scales.L_turb1.*sqrt(tke1),zunitsp,'b',...
         work.kdifft1,zunitsp,'k',...
         work.kdifft1x,zunitsp,'k--')
    ylim([0,zplottop])
    title('Diffusivity')
    xlabel('K_1')
    ylabel(labelz)
    set(gca,'FontSize',fs)
else
    figure(5)
    subplot(1,1,1)
    plot(scales.L_turb1.*sqrt(tke1),zunitsp,'b',...
         work.kdifft1,zunitsp,'k')
    ylim([0,zplottop])
    title('Diffusivity')
    xlabel('K_1')
    ylabel(labelz)
    set(gca,'FontSize',fs)
    pause
end
if plottype == 0
    figure(4)
    subplot(2,4,8)
    plot(scales.L_turb2.*sqrt(tke2),zunitsp,'b',...
         work.kdifft2,zunitsp,'k',...
         work.kdifft2x,zunitsp,'k--')
    ylim([0,zplottop])
    title('Diffusivity')
    xlabel('K_2')
    ylabel(labelz)
    set(gca,'FontSize',fs)
else
    figure(5)
    subplot(1,1,1)
    plot(scales.L_turb2.*sqrt(tke2),zunitsp,'b',...
         work.kdifft2,zunitsp,'k')
    ylim([0,zplottop])
    title('Diffusivity')
    xlabel('K_2')
    ylabel(labelz)
    set(gca,'FontSize',fs)
    pause
end

