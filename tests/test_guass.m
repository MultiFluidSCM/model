% Test Gaussian cloud model formulas

rr2pi = 1/sqrt(2*pi);
rr2 = 1/sqrt(2);

N = 1000;
xsmax = 5;
xs = linspace(-xsmax,xsmax,N+1);
dx = 2*xsmax/N;

for k = 1:N+1
    
    % Brute force integrals
    R = 0;
    ql = 0;
    qlql =0;
    for j = 1:k-1
        G1 = rr2pi*exp(-xs(j  )*xs(j  )/2);
        G2 = rr2pi*exp(-xs(j+1)*xs(j+1)/2);
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
    ra(k) = 0.5*(1 + erf(-xs(k)*rr2));
    I1 = rr2pi*exp(-xs(k)*xs(k)/2));
    qla(k) = -xs(k)*ra(k) + I1;
    qlqla(k) = ra(k)*(1 + xs(k)^2) ...
             + (ra(k) - 1)*xs(k)*I1 ...
             - ra(k)^2*xs(k)^2 ...
             - I1*I1;
         
end

% Plot

figure(1)
subplot(3,1,1)
plot(xs,rb,'b',xs,ra,'r')
set(gca,'fontsize',18)
xlabel('xs')
ylabel('R')

subplot(3,1,2)
plot(xs,qlb,'b',xs,qla,'r')
set(gca,'fontsize',18)
xlabel('xs')
ylabel('ql')

subplot(3,1,1)
plot(xs,rb,'b',xs,ra,'r')
set(gca,'fontsize',18)
xlabel('xs')
ylabel('Var ql')    

    
end


