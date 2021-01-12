% Test Gaussian cloud model formulas


N = 1000;
xsmax = 5;
xs = linspace(-xsmax,xsmax,N+1);
dx = 2*xsmax/N;


for pass = 1:2

for k = 1:N+1
    
    % Brute force integrals
    R = 0;
    ql = 0;
    qlql =0;
    for j = k:N
        if pass == 1
            G1 = G_gauss(xs(j  ));
            G2 = G_gauss(xs(j+1));
        else
            G1 = G_uniform(xs(j  ));
            G2 = G_uniform(xs(j+1));
        end
        R = R + 0.5*(G1 + G2)*dx;
        ql = ql + 0.5*((xs(j  ) - xs(k))*G1 ...
                     + (xs(j+1) - xs(k))*G2)*dx;
        qlql = qlql + 0.5*((xs(j  ) - xs(k))^2*G1 ...
                         + (xs(j+1) - xs(k))^2*G2)*dx;
    end
    qlql = qlql - ql*ql;
    rb(k) = R;
    qlb(k) = ql;
    qlqlb(k) = qlql;
    
    % Formulas
    if pass == 1
        [I0, I1, I2] = I_gauss(xs(k));
    else
        [I0, I1, I2] = I_uniform(xs(k));
    end
    ra(k) = I0;
    qla(k) = -xs(k)*I0 + I1;
%     qlqla(k) = ra(k)*(1 + xs(k)^2) ...
%              + (2*ra(k) - 1)*xs(k)*I1 ...
%              - ra(k)^2*xs(k)^2 ...
%              - I1*I1;
%              - qla(k)*qla(k);
    qlqla(k) = I2 - 2*xs(k)*I1 + xs(k)*xs(k)*I0 - qla(k)^2;
end

% Plot

xplotmax = 3;

figure(1)
subplot(4,1,1)
plot(xs,rb,'b+',xs,ra,'r')
xlim([-xplotmax,xplotmax])
set(gca,'fontsize',18)
xlabel('xs')
ylabel('R')
if pass == 1
    hold on
else
    hold off
end

subplot(4,1,2)
plot(xs,qlb,'b+',xs,qla,'r')
xlim([-xplotmax,xplotmax])
set(gca,'fontsize',18)
xlabel('xs')
ylabel('ql')
if pass == 1
    hold on
else
    hold off
end

subplot(4,1,3)
plot(xs,qlqlb,'b+',xs,qlqla,'r')
xlim([-xplotmax,xplotmax])
set(gca,'fontsize',18)
xlabel('xs')
ylabel('Var ql')
if pass == 1
    hold on
else
    hold off
end

subplot(4,1,4)
plot(xs,qlqla./(qla.*qla),'k')
xlim([-xplotmax,xplotmax])
ylim([0,15])
set(gca,'fontsize',18)
xlabel('xs')
ylabel('qlql/ql^2')
if pass == 1
    hold on
else
    hold off
end

end