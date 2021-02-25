% Plot time series of maximum residuals as a check on convergence

% Figure 27

figure(27)

if exist('qn_iter_max')
    ts.maxres1w(plotstep) = maxres1w(qn_iter_max);
    ts.maxres2w(plotstep) = maxres2w(qn_iter_max);
    ts.maxres1m(plotstep) = maxres1m(qn_iter_max);
    ts.maxres2m(plotstep) = maxres2m(qn_iter_max);
    ts.maxres1eta(plotstep) = maxres1eta(qn_iter_max);
    ts.maxres2eta(plotstep) = maxres2eta(qn_iter_max);
    ts.maxres1q(plotstep) = maxres1q(qn_iter_max);
    ts.maxres2q(plotstep) = maxres2q(qn_iter_max);
    ts.maxres1u(plotstep) = maxres1u(qn_iter_max);
    ts.maxres2u(plotstep) = maxres2u(qn_iter_max);
    ts.maxres1v(plotstep) = maxres1v(qn_iter_max);
    ts.maxres2v(plotstep) = maxres2v(qn_iter_max);
    ts.maxres1eosetap(plotstep) = maxres1eosetap(qn_iter_max);
    ts.maxres2eosetap(plotstep) = maxres2eosetap(qn_iter_max);
    ts.maxres1eoseta(plotstep) = maxres1eoseta(qn_iter_max);
    ts.maxres2eoseta(plotstep) = maxres2eoseta(qn_iter_max);
    ts.maxressigma(plotstep) = maxressigma(qn_iter_max);
    ts.maxres1tke(plotstep) = maxres1tke(qn_iter_max);
    ts.maxres2tke(plotstep) = maxres2tke(qn_iter_max);
    ts.maxres1Veta(plotstep) = maxres1Veta(qn_iter_max);
    ts.maxres2Veta(plotstep) = maxres2Veta(qn_iter_max);
    ts.maxres1Vq(plotstep) = maxres1Vq(qn_iter_max);
    ts.maxres2Vq(plotstep) = maxres2Vq(qn_iter_max);
else
    ts.maxres1w(plotstep) = 0;
    ts.maxres2w(plotstep) = 0;
    ts.maxres1m(plotstep) = 0;
    ts.maxres2m(plotstep) = 0;
    ts.maxres1eta(plotstep) = 0;
    ts.maxres2eta(plotstep) = 0;
    ts.maxres1q(plotstep) = 0;
    ts.maxres2q(plotstep) = 0;
    ts.maxres1u(plotstep) = 0;
    ts.maxres2u(plotstep) = 0;
    ts.maxres1v(plotstep) = 0;
    ts.maxres2v(plotstep) = 0;
    ts.maxres1eosetap(plotstep) = 0;
    ts.maxres2eosetap(plotstep) = 0;
    ts.maxres1eoseta(plotstep) = 0;
    ts.maxres2eoseta(plotstep) = 0;
    ts.maxressigma(plotstep) = 0;
    ts.maxres1tke(plotstep) = 0;
    ts.maxres2tke(plotstep) = 0;
    ts.maxres1Veta(plotstep) = 0;
    ts.maxres2Veta(plotstep) = 0;
    ts.maxres1Vq(plotstep) = 0;
    ts.maxres2Vq(plotstep) = 0;
end

subplot(3,3,1)
plot(ts.time/3600,ts.maxres1w,'b',ts.time/3600,ts.maxres2w,'r')
title('Res w')
xlim([0,time.tstop/3600])
set(gca,'fontsize',fs)

subplot(3,3,2)
plot(ts.time/3600,ts.maxres1m,'b',ts.time/3600,ts.maxres2m,'r')
title('Res m')
xlim([0,time.tstop/3600])
set(gca,'fontsize',fs)

subplot(3,3,3)
plot(ts.time/3600,ts.maxres1eta,'b',ts.time/3600,ts.maxres2eta,'r')
title('Res eta')
xlim([0,time.tstop/3600])
set(gca,'fontsize',fs)

subplot(3,3,4)
plot(ts.time/3600,ts.maxres1q,'b',ts.time/3600,ts.maxres2q,'r')
title('Res q')
xlim([0,time.tstop/3600])
set(gca,'fontsize',fs)

subplot(3,3,5)
plot(ts.time/3600,ts.maxres1u,'b',ts.time/3600,ts.maxres2u,'r',...
     ts.time/3600,ts.maxres1v,'b',ts.time/3600,ts.maxres2v,'r')
title('Res u, v')
xlim([0,time.tstop/3600])
set(gca,'fontsize',fs)

subplot(3,3,6)
plot(ts.time/3600,ts.maxres1eosetap,'b',ts.time/3600,ts.maxres2eosetap,'r',...
     ts.time/3600,ts.maxres1eoseta,'b--',ts.time/3600,ts.maxres2eoseta,'r--',...
     ts.time/3600,1e4*ts.maxressigma,'k')
title('Res eos')
xlim([0,time.tstop/3600])
set(gca,'fontsize',fs)

subplot(3,3,7)
plot(ts.time/3600,ts.maxres1tke,'b',ts.time/3600,ts.maxres2tke,'r')
title('Res tke')
xlim([0,time.tstop/3600])
set(gca,'fontsize',fs)

subplot(3,3,8)
plot(ts.time/3600,ts.maxres1Veta,'b',ts.time/3600,ts.maxres2Veta,'r')
title('Res Veta')
xlim([0,time.tstop/3600])
set(gca,'fontsize',fs)

subplot(3,3,9)
plot(ts.time/3600,ts.maxres1Vq,'b',ts.time/3600,ts.maxres2Vq,'r')
title('Res Vq')
xlim([0,time.tstop/3600])
set(gca,'fontsize',fs)




