% Check the linearization of the sorting detrainment factor chi_hat

clear

rr2 = sqrt(0.5);
r2bypi = sqrt(2/pi);


w2 = 1;
sigmaw = [0.01:0.01:5];
% chi_cut = max(-rr2*w2./sigmaw,-100);
% E = exp(-chi_cut.^2);
% R = (1 + erf(chi_cut));
% alt = chi_cut < -4;
% chi_hat = (1 - alt).*(-r2bypi*E./(R + alt)) ...
%         + alt.*sqrt(2).*chi_cut./(1 - 0.5./chi_cut.^2);
for k = 1:numel(sigmaw)
    
    chi_cut(k) = -rr2*w2/sigmaw(k);
    [a,b] = compute_chi_hat(chi_cut(k));
    chi_hat(k) = a;
    formula(k) = b;
    
end
what12 = w2 + sigmaw.*chi_hat;

subplot(3,2,1)
% plot(sigmaw,chi_hat,'b',sigmaw,sigmaw.*chi_hat,'r',sigmaw,sqrt(2)*chi_cut,'g')
plot(sigmaw,sigmaw.*chi_hat,'r',sigmaw,sqrt(2)*sigmaw.*chi_cut,'g')
title('chi-hat and sigma*chi-hat')

subplot(3,2,3)
plot(sigmaw,what12,'b')
title('what12')

dwhatdsigma = (what12(2:end) - what12(1:end-1))./(sigmaw(2:end) - sigmaw(1:end-1));
s = 0.5*(sigmaw(2:end) + sigmaw(1:end-1));
%formula = (-r2bypi*R.*E - 2*r2bypi*(chi_cut.^2).*R.*E - (2*sqrt(2)/pi)*chi_cut.*E.*E)./(R.*R);
%formula = (1 - alt).*formula + alt.*0;
subplot(3,2,5)
plot(s,dwhatdsigma,'r',sigmaw,formula,'ko')

w2 = -1;
sigmaw = [0.01:0.01:5];
% chi_cut = max(-rr2*w2./sigmaw,-100);
% E = exp(-chi_cut.^2);
% R = (1 + erf(chi_cut));
% alt = chi_cut < -4;
% chi_hat = (1 - alt).*(-r2bypi*E./(R + alt)) ...
%         + alt.*sqrt(2).*chi_cut./(1 - 0.5./chi_cut.^2);
for k = 1:numel(sigmaw)
    
    chi_cut(k) = -rr2*w2/sigmaw(k);
    [a,b] = compute_chi_hat(chi_cut(k));
    chi_hat(k) = a;
    formula(k) = b;
    
end
what12 = w2 + sigmaw.*chi_hat;

subplot(3,2,2)
plot(sigmaw,chi_hat,'b',sigmaw,sigmaw.*chi_hat,'r')
title('chi-hat and sigma*chi-hat')

subplot(3,2,4)
plot(sigmaw,what12,'b')
title('what12')

dwhatdsigma = (what12(2:end) - what12(1:end-1))./(sigmaw(2:end) - sigmaw(1:end-1));
s = 0.5*(sigmaw(2:end) + sigmaw(1:end-1));
%formula = (-r2bypi*R.*E - 2*r2bypi*(chi_cut.^2).*R.*E - (2*sqrt(2)/pi)*chi_cut.*E.*E)./(R.*R);
subplot(3,2,6)
plot(s,dwhatdsigma,'r',sigmaw,formula,'ko')