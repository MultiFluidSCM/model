function [Vartheta] = theta_variance(eta,q,vareta,varq,Varql,Covarqleta,Covarqlq,therm,p,p00)
%THETA_VARIANCE Compute potential temperature variance
%   
% Newton iteration to compute liquid water potential temperature thetal
% It's just the Gibbs function derivatives we need
tt = 300.0;
for iter = 1:4
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a ] = gibbs(p00,tt,q,therm);
    res = eta + gt;
    t_inc = -res/gtt;
    tt = tt + t_inc;
end
thetal = tt;

% Exner function
exner = (p/p00)^therm.kappa;

% Covariance of eta and q
rr = 0;
covaretaq = rr*sqrt(vareta*varq);

% Theta variance
Vartheta = vareta/(gtt^2) ...
         + Varql*(therm.L0/(exner*therm.Cpd))^2 ...
         + varq *(gtw/gtt)^2 ...
         - Covarqleta*(2*therm.L0/(exner*therm.Cpd*gtt)) ...
         - Covarqlq  *(2*therm.L0*gtw/(exner*therm.Cpd*gtt)) ...
         + covaretaq *(2*gtw/gtt^2);

end

