% Estimate and plot profiles of standard deviation of various
% quantities

figure(20)



% Estimate standard deviation of tracers as a combination of two terms:
% (i) Turbulent length scale * vertical derivative
% (ii) Difference between fluids 1 and 2
% deta2dz = (eta2(2:nzp) - eta2(1:nz))./dzp;
% dq2dz   = (q2(2:nzp)   - q2(1:nz))  ./dzp;
% term1 = scales.L_turb2.*deta2dz;
% term2 = eta2bar - eta1bar;
% etastd0 = sqrt(term1.^2 + term2.^2);
% term1 = scales.L_turb2.*dq2dz;
% term2 = q2bar - q1bar;
% qstd0 = sqrt(term1.^2 + term2.^2);

% Values from more complete calculation
etastd1 = weight_to_w(grid,sqrt(state_new.fluid(1).vareta));
qstd1   = weight_to_w(grid,sqrt(state_new.fluid(1).varq  ));
etastd2 = weight_to_w(grid,sqrt(state_new.fluid(2).vareta));
qstd2   = weight_to_w(grid,sqrt(state_new.fluid(2).varq  ));
qlstd1  = sqrt(max(0,eos.Varql1));
qlstd2  = sqrt(max(0,eos.Varql2));
bstd    = constants.phys.gravity*sqrt(eos.Varrhow2)./eos.rhow2;

% Subfilter standard deviation of w, assuming tke has equal contributions
% from u, v and w
wstd1 = sqrt(tke1*2/3);
wstd2 = sqrt(tke2*2/3);

% Possible normalizations for b for candidate smooth switches
% bnorm = bstd;
temp = tke2./scales.L_turb2;
bnorm = weight_to_w(grid,temp);

% And plot
subplot(3,4,1)
plot(wstd1,zunitsp,'b',wstd2,zunitsp,'r')
set(gca,'FontSize',fs)
ylim([0,zplottop])
title('std dev w ')
xlabel('\sigma_w')
ylabel(labelz)

subplot(3,4,2)
plot(w2bar - wstd2,zunitsp,'r',...
     w2bar + wstd2,zunitsp,'r',...
     state_new.fluid(1).w,zunitsw,'b')
set(gca,'FontSize',fs)
ylim([0,zplottop])
title('w_2 spread')
xlabel('w')
ylabel(labelz)

subplot(3,4,3)
plot(etastd1,zunitsw,'b',etastd2,zunitsw,'r')
set(gca,'FontSize',fs)
ylim([0,zplottop])
title('std dev eta ')
xlabel('\sigma_\eta')
ylabel(labelz)

subplot(3,4,4)
plot(eta2 - etastd2,zunitsw,'r',...
     eta2 + etastd2,zunitsw,'r',...
     state_new.fluid(1).eta,zunitsw,'b')
set(gca,'FontSize',fs)
ylim([0,zplottop])
title('eta_2 spread')
xlabel('\eta')
ylabel(labelz)

subplot(3,4,5)
plot(qstd1,zunitsw,'b',qstd2,zunitsw,'r')
set(gca,'FontSize',fs)
ylim([0,zplottop])
title('std dev q ')
xlabel('\sigma_q')
ylabel(labelz)

subplot(3,4,6)
plot(q2 - qstd2,zunitsw,'r',...
     q2 + qstd2,zunitsw,'r',...
     state_new.fluid(1).q,zunitsw,'b')
set(gca,'FontSize',fs)
ylim([0,zplottop])
title('q_2 spread')
xlabel('q')
ylabel(labelz)

subplot(3,4,7)
plot(qlstd1,zunitsw,'b',...
     qlstd2,zunitsw,'r')
ylim([0,zplottop])
title('std dev ql')
xlabel('\sigma_{ql}')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(3,4,8)
plot(eos.ql2 - qlstd2,zunitsw,'r--',...
     eos.ql2 + qlstd2,zunitsw,'r--',...
     eos.ql2,zunitsw,'r')
set(gca,'FontSize',fs)
ylim([0,zplottop])
title('ql_2 spread')
xlabel('ql')
ylabel(labelz)

subplot(3,4,9)
plot(bstd,zunitsw,'k',bnorm,zunitsw,'r')
ylim([0,zplottop])
title('std dev buoy, b-norm')
xlabel('\sigma_b')
ylabel(labelz)
set(gca,'FontSize',fs)

subplot(3,4,10)
plot(work.buoy - bstd,zunitsw,'k--',...
     work.buoy + bstd,zunitsw,'k--',...
     work.buoy,zunitsw,'k')
set(gca,'FontSize',fs)
ylim([0,zplottop])
title('b spread')
xlabel('b')
ylabel(labelz)

% % Candidate smooth switches
% ss1 = 0.5*(1 - tanh(   work.buoy./bnorm + 2));
% ss2 = 0.5*(1 - tanh( 2*work.buoy./bnorm + 2));
% ss3 = 0.5*(1 - tanh( 5*work.buoy./bnorm + 2));
% ss4 = 0.5*(1 - tanh(10*work.buoy./bnorm + 2));
% 
% subplot(3,4,11)
% plot(ss1,zunitsw,'k--',...
%      ss2,zunitsw,'b',...
%      ss3,zunitsw,'k',...
%      ss4,zunitsw,'r',...
%      relabel.ss3bar,zunitsp,'o')
% set(gca,'FontSize',fs)
% ylim([0,zplottop])
% title('Smooth switches')
% xlabel('Smooth switches')
% ylabel(labelz)

subplot(3,4,11)
plot(state_new.fluid(1).vareta.*(eos.theta1p/constants.therm.Cpd).^2,zunitsp,'k',...
     eos.Vartheta1,zunitsw,'r',...
     eos.Vartheta1_alt,zunitsw,'b')
set(gca,'FontSize',fs)
ylim([0,zplottop])
title('Theta1 variance')
xlabel('Theta variance')
ylabel(labelz)

subplot(3,4,12)
plot(state_new.fluid(2).vareta.*(eos.theta2p/constants.therm.Cpd).^2,zunitsp,'k',...
     eos.Vartheta2,zunitsw,'r',...
     eos.Vartheta2_alt,zunitsw,'b')
set(gca,'FontSize',fs)
ylim([0,zplottop])
title('Theta2 variance')
xlabel('Theta variance')
ylabel(labelz)
