% Plot the cubic fitted to find zstar using option 4

figure(18)
hold off
kk1 = kmxdtdz-1;
kk2 = kmxdtdz+2;
kk1 = 1;
kk2 = 65; %10;
plot(Tw(kk1:kk2),zw(kk1:kk2),'k-o')
hold on

nplt = 51;
zxx = linspace(zw(kmxdtdz-1),zw(kmxdtdz+2),nplt);
a = dz*(dTp + dT0) - 2*(Tp - T0);
b = 3*(Tp - T0) - dz*(dTp + 2*dT0);
c = dz*dT0;
d = T0;
Tstar = a*zhat^3 + b*zhat^2 + c*zhat + d;
plot([Tstar],[zstar4],'r*')

for k = 1:nplt
    zhat = (zxx(k) - zz0)/dz;
    Txx(k) = a*zhat^3 + b*zhat^2 + c*zhat + d;
end
plot(Txx,zxx,'b')

hold off
xlabel('T')
ylabel('z')
%xlim([296.8,298.6])
pause

figure(1)
