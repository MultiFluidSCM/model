% Plot time series of surface forcing for ARM case

% Set information needed for forcing
% Set time (x) and sensible (shf) - latent heat fluxes (lhf) in W m-2
x = [0; 14400; 23400; 27000; 36000; 45000; 52200];
shf2 = [1; 90; 140; 140; 100; 1; -10];
shf = [-30; 90; 140; 140; 100; -10; -10];
lhf = [5; 250; 450; 500; 420; 180; 0];
       
x = x/3600;       

figure(6)
subplot(3,1,1)
fs = 12;
hold off
plot(x,shf,'b','linewidth',1.5)
hold on
plot(x,shf2,'b--','linewidth',1.5)
plot(x,lhf,'r','linewidth',1.5)
title('ARM sfc forcing','fontsize',12)
xlabel('Time (hrs)','fontsize',12)
ylabel('SHF + LHF (W/m^2)','Fontsize',12)
xlim([0,x(end)])

