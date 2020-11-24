% Plot profiles of fluid 2 cloud fraction

% Assume cloud fractions have been computed already in find_eos_sg

figure(21)

% Cloud frcation and updraft fraction
subplot(1,3,1)
plot(cld_frac1,zunitsw,'b',...
     cld_frac2,zunitsw,'r',...
     tot_cld_frac,zunitsw,'k',...
     sigma2,zunitsp,'k--')
set(gca,'FontSize',fs)
%axis([0,0.4,0,zplottop])
ylim([0,zplottop])
title('Cloud frac')
xlabel('f')
ylabel(labelz)

% Cloud fraction only
subplot(1,3,2)
plot(cld_frac1,zunitsw,'b',...
     cld_frac2,zunitsw,'r',...
     tot_cld_frac,zunitsw,'k')
set(gca,'FontSize',fs)
%axis([0,0.4,0,zplottop])
ylim([0,zplottop])
title('Cloud frac')
xlabel('f')
ylabel(labelz)
