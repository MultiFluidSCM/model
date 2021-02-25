% Compare turbulence inverse timescale with sorting rate N

figure(22)

n1pos = sqrt(max( eos.nsq1,0));
dw2dz = (w2(2:nzp) - w2(1:nz))./dzp;
rTbuoy = -(abovep.*work.buoy(2:nzp) + belowp.*work.buoy(1:nz))./sqrt(tke2);

plot(1./scales.T_turb2,zunitsp,'k',n1pos,zunitsp,'b',max(-dw2dz,0),zunitsp,'r',rTbuoy,zunitsp,'b--')
set(gca,'FontSize',fs)
legend('r_{turb2}','N1^+','-dw2dz^+','r_{buoy}')
ylim([0,zplottop])
xlim([0,0.05])
title(' Rates ')
xlabel(' 1/T, N, -dw2/dz, buoy-rate ')
ylabel(labelz)

