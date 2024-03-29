% Plot basic fields

% Figure 1: Basic fields mass frac, w, eta, q, u and v, buoyancy, N^2,
% thetal, RH, tke


% Updraft mass fraction and target
massfrac = m2./(m1 + m2);
if plottype == 0
    figure(1)
    subplot(2,5,1)
    plot(massfrac,zunitsp,'k',relabel.ideal,zunitsp,'k--')
    %axis([0,0.25,0,zplottop])
    ylim([0,zplottop])
    title('Updraft fraction')
    xlabel('Updraft frac')
    ylabel(labelz)
    set(gca,'FontSize',fs)
else
    figure(5)
    subplot(1,1,1)
    plot(massfrac,zunitsp,'k','Linewidth',1.5)
%    axis([0,0.25,0,zplottop])
    axis([0,1,0,zplottop])
    title('(a) Updraft frac')
    xlabel('\sigma_2')
    ylabel(labelz)
    set(gca,'FontSize',fs)
    pause
end


% Updraft and environment vertical velocity
if plottype == 0
    figure(1)
    subplot(2,5,2)
%     plot(w2/wstar,zunitsw,'k',w1/wstar,zunitsw,'k--')
%     axis([-1,1.5,0,zplottop])
%     title(['w_i   w_* = ',num2str(wstar,'%0.3g')])
%     xlabel('w_i / w_*')
%     ylabel(labelz)
    plot(w2,zunitsw,'k',w1,zunitsw,'k--')
    axis([-1,2.0,0,zplottop])
    title(['w_i   w_* = ',num2str(wstar,'%0.3g')])
    xlabel('w_i')
    ylabel(labelz)
    set(gca,'FontSize',fs)
else
    figure(5)
%     plot(w2/wstar,zunitsw,'k',w1/wstar,zunitsw,'k--')
%     axis([-1,1.5,0,zplottop])
%     title('(b) Scaled updraft w ')
%     xlabel('w_2 / w_*')
%     ylabel(labelz)
    plot(w2,zunitsw,'k',w1,zunitsw,'k--')
    axis([-1,2.0,0,zplottop])
    % title('(b) Updraft w ')
    title('Up- and downdraft w')
    xlabel('w_i')
    ylabel(labelz)
    set(gca,'FontSize',fs)
    pause
end


% Updraft and environment eta
if plottype == 0
    figure(1)
    subplot(2,5,3)
    plot(eta2,zunitsw,'k',eta1,zunitsw,'k--')
    ylim([0,zplottop])
    title('\eta_1 and \eta_2 ')
    xlabel('\eta_1 and \eta_2')
    ylabel(labelz)
    set(gca,'FontSize',fs)
else
    figure(5)
    plot(eta2,zunitsw,'k','LineWidth',1.5)
    hold on
    plot(eta1,zunitsw,'k--','LineWidth',1.5)
    hold off
    ylim([0,zplottop])
    title('\eta_1 and \eta_2 ')
    xlabel('\eta_1 and \eta_2')
    ylabel(labelz)
    set(gca,'FontSize',fs)
    pause
end


% Liquid water and vapour profiles
% and relative humidity
for k = 1:nzp
    if k == 1
        pbar = grid.extrapb1*p(1) + grid.extrapb2*p(2);
    elseif k == nzp
        pbar = grid.extraptnz*p(nz) + grid.extraptnzm*p(nz-1);
    else
        pbar   = grid.abovew(k)*p(k) ...
               + grid.beloww(k)*p(k-1);
    end
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(pbar,Tw1(k),q1(k),constants.therm);
    vapour1(k) = (1 - q1(k))*(1 - a)/a;
    liquid1(k) = q1(k) - vapour1(k);
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(pbar,Tw2(k),q2(k),constants.therm);
    vapour2(k) = (1 - q2(k))*(1 - a)/a;
    liquid2(k) = q2(k) - vapour2(k);
    rh1(k) = relhum(pbar,Tw1(k),q1(k),constants.therm);
    rh2(k) = relhum(pbar,Tw2(k),q2(k),constants.therm);
end
disp(['Max liquid ',num2str(max(liquid2))])
if plottype == 0
    figure(1)
    subplot(2,5,4)
    plot(vapour1,zunitsw,'b--',...
         100*liquid1,zunitsw,'r--',...
         vapour2,zunitsw,'b',...
         100*liquid2,zunitsw,'r',...  %)
         100*eos.ql1,zunitsw,'g',...
         100*eos.ql2,zunitsw,'k')
    ylim([0,zplottop])
    %xlim([0,20e-3])
    title('Vapour + liquid ')
    xlabel('q, 100q_l')
    ylabel(labelz)
    set(gca,'FontSize',fs)
else
    figure(5)
    %if ipt == 1
    %        hold off
    %end
    plot(vapour1,zunitsw,'b--','LineWidth',1.5)
    hold on
    plot(100*liquid1,zunitsw,'r--','LineWidth',1.5)
    plot(vapour2,zunitsw,'b','LineWidth',1.5)
    plot(100*liquid2,zunitsw,'r','LineWidth',1.5)
    hold off
    ylim([0,zplottop])
    xlim([0,20e-3])
    title('Vapour + liquid ')
    xlabel('q, 100q_l')
    ylabel(labelz)
    set(gca,'FontSize',fs)
    pause
end

rhmax = round(100*max(rh2));
if plottype == 0
    subplot(2,5,9)
    plot(rh1,zunitsw,'b',rh2,zunitsw,'r')
    axis([0,1.1,0,zplottop])
    title(['RH ',num2str(rhmax),'%'])
    xlabel('RH')
    ylabel(labelz)
    set(gca,'FontSize',fs)
else
    figure(5)
    plot(rh1,zunitsw,'b','Linewidth',1.5)
    plot(rh2,zunitsw,'r','Linewidth',1.5)
    axis([0,1.1,0,zplottop])
    title('Rel. Hum. ')
    xlabel('RH')
    ylabel(labelz)
    set(gca,'FontSize',fs)
    pause
end


% u and v profiles
if plottype == 0
    subplot(2,5,5)
    plot(u1,zunitsp,'b',u2,zunitsp,'r',v1,zunitsp,'b--',v2,zunitsp,'r--')
    axis([-12,12,0,zplottop])
    title('u and v')
    xlabel('u and v')
    ylabel(labelz)
    set(gca,'FontSize',fs);
else
    figure(5)
    plot(u1,zunitsp,'b',u2,zunitsp,'r',v1,zunitsp,'b--',v2,zunitsp,'r--')
    axis([-12,12,0,zplottop])
    title('u and v')
    xlabel('u and v')
    ylabel(labelz)
    set(gca,'FontSize',fs);
    pause
end


% Buoyancy of updraft
if plottype == 0
    figure(1)
    subplot(2,5,6)
    plot(work.buoy,zunitsw,'k')
    ylim([0,zplottop])
    title('Updraft buoyancy ')
    xlabel('b')
    ylabel(labelz)
    set(gca,'FontSize',fs)
else
    figure(5)
    plot(work.buoy,zunitsw,'k','linewidth',1.5)
    ylim([0,zplottop])
    title('Updraft buoyancy ')
    xlabel('b')
    ylabel(labelz)
    set(gca,'FontSize',fs)
    pause
end


% Buoyancy frequency
n1pos = sqrt(max( eos.nsq1,0));
n1neg = sqrt(max(-eos.nsq1,0));
%n2posx = sqrt(max( work.nsq2,0));
%n2negx = sqrt(max(-work.nsq2,0));
n2pos = sqrt(max( eos.nsq2,0));
n2neg = sqrt(max(-eos.nsq2,0));
if plottype == 0
    subplot(2,5,7)
    plot(n1pos,zunitsp,'b--',n1neg,zunitsp,'r--',n2pos,zunitsp,'b',n2neg,zunitsp,'r')
    %plot(n1pos,zunitsp,'b--',n1neg,zunitsp,'r--',n2pos,zunitsp,'b',n2neg,zunitsp,'r',...
    %     n2posx,zunitsp,'kx',n2negx,zunitsp,'ko')
%    plot(n1pos,zunitsp,'b',n1neg,zunitsp,'r',n1posx,zunitsp,'kx',n1negx,zunitsp,'ko')
    axis([0,0.04,0,zplottop])
    title('N or (-N^2)^{1/2}')
    xlabel('N')
    ylabel(labelz)
    set(gca,'FontSize',fs);
else
    figure(5)
    plot(n1pos,zunitsp,'b--',n1neg,zunitsp,'r--',n2pos,zunitsp,'b',n2neg,zunitsp,'r')
    axis([0,0.04,0,zplottop])
    title('N or (-N^2)^{1/2}')
    xlabel('N')
    ylabel(labelz)
    set(gca,'FontSize',fs);
    pause
end


% Updraft and environment (density) potential temperature
if plottype == 0
    figure(1)
    subplot(2,5,8)
    plot(eos.theta2,zunitsw,'k',eos.theta1,zunitsw,'k--')
    % plot(eos.theta_rho2,zunitsw,'k',eos.theta_rho1,zunitsw,'k--')
    ylim([0,zplottop])
    %title('\theta_{\rho 1} and \theta_{\rho 2} ')
    %xlabel('\theta_{\rho 1} and \theta_{\rho 2}')
    title('\theta_{1} and \theta_{2} ')
    xlabel('\theta_{1} and \theta_{2}')
    ylabel(labelz)
    set(gca,'FontSize',fs)
else
    figure(5)
    plot(eos.theta2,zunitsw,'k','LineWidth',1.5)
    hold on
    plot(eos.theta1,zunitsw,'k--','LineWidth',1.5)
    hold off
%     plot(eos.theta_rho2,zunitsw,'k','LineWidth',1.5)
%     hold on
%     plot(eos.theta_rho1,zunitsw,'k--','LineWidth',1.5)
%     hold off
    ylim([0,zplottop])
    %title('\theta_{\rho 1} and \theta_{\rho 2} ')
    %xlabel('\theta_{\rho 1} and \theta_{\rho 2}')
    title('\theta_{1} and \theta_{2} ')
    xlabel('\theta_{1} and \theta_{2}')
    ylabel(labelz)
    set(gca,'FontSize',fs)
    pause
end


% % Updraft and environment potential density
% if plottype == 0
%     figure(1)
%     subplot(2,5,8)
%     plot(eos.rho_ptl2,zunitsw,'k',eos.rho_ptl1,zunitsw,'k--')
%     ylim([0,zplottop])
%     title('\rho_{1 ptl} and \rho_{2 ptl}')
%     xlabel('\rho_{1 ptl} and \rho_{2 ptl}')
%     ylabel(labelz)
%     set(gca,'FontSize',fs)
% else
%     figure(5)
%     plot(eos.rho_ptl2,zunitsw,'k','LineWidth',1.5)
%     hold on
%     plot(eos.rho_ptl1,zunitsw,'k--','LineWidth',1.5)
%     hold off
%     ylim([0,zplottop])
%     title('\rho_{1 ptl} and \rho_{2 ptl}')
%     xlabel('\rho_{1 ptl} and \rho_{2 ptl}')
%     ylabel(labelz)
%     set(gca,'FontSize',fs)
%     pause
% end


% TKE
% normalization
tkenorm = wstar*wstar;
% KE of mean resolved
umean = (m1.*u1 + m2.*u2)./(m1 + m2);
vmean = (m1.*v1 + m2.*v2)./(m1 + m2);
wmean = (m1bar.*w1 + m2bar.*w2)./(m1bar + m2bar);
wpert1sq = (wmean - w1).^2;
wpert2sq = (wmean - w2).^2;
% Resolved TKE per unit mass
Rtke1 = 0.5*((umean - u1).^2 + (vmean - v1).^2 + grid.aboves.*wpert1sq(2:nzp) + grid.belows.*wpert1sq(1:nz));
Rtke2 = 0.5*((umean - u2).^2 + (vmean - v2).^2 + grid.aboves.*wpert2sq(2:nzp) + grid.belows.*wpert2sq(1:nz));
% Total TKE per unit mass
tketot = (m1.*(Rtke1 + tke1) + m2.*(Rtke2 + tke2))./(m1 + m2);
if plottype == 0
    subplot(2,5,10)
    plot(tke1 /tkenorm,zunitsp,'b'  ,tke2 /tkenorm,zunitsp,'r',  ...
         Rtke1/tkenorm,zunitsp,'b--',Rtke2/tkenorm,zunitsp,'r--', ...
         tketot/tkenorm,zunitsp,'k')
    ylim([0,zplottop])
    title('TKE 1, 2 and tot')
    xlabel('TKE/w_*^2 ')
    ylabel(labelz)
    set(gca,'FontSize',fs)
else
    figure(5)
    plot(tke1 /tkenorm,zunitsp,'b'  ,tke2 /tkenorm,zunitsp,'r',  ...
         Rtke1/tkenorm,zunitsp,'b--',Rtke2/tkenorm,zunitsp,'r--', ...
         tketot/tkenorm,zunitsp,'k')
    ylim([0,zplottop])
    title('TKE 1, 2 and tot')
    xlabel('TKE/w_*^2 ')
    ylabel(labelz)
    set(gca,'FontSize',fs)
    pause
end

