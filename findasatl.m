function [ asat ] = findasatl( p, T, ww, therm)
%FINDASATL Compute dry mass fraction of saturated air
%   Determine the mass fraction of dry air in the gaseous
%   part of a moist air - liquid water mix, given
%   pressure, temperature, and mass fraction of total water.
%   A Newton method is used.
%   

% Unpack constants needed locally
p0v = therm.p0v;
T0  = therm.T0;
Cpv = therm.Cpv;
Cl  = therm.Cl;
Rv  = therm.Rv;
L00 = therm.L00;
epsilon = therm.epsilon;

% First guess
psat = p0v*(T/T0)^((Cpv-Cl)/Rv)*exp((L00/Rv)*(1.0d0/T0 - 1.0d0/T));
asat = (p - psat)/(p + (epsilon - 1.0d0)*psat);


% Chemical potential (=Gibbs) of the liquid; only need to
% compute this once
[mul,gp,gt,gpp,gpt,gtt] = gibbsl(p,T,therm);

for iter = 1:3
    % Chemical potential of the vapour in the air-vapour mix
    [g,gp,gt,ga,gpp,gpt,gpa,gtt,gta,gaa,gtw] = gibbsav(p,T,asat,therm);
    muav = g - asat*ga;
    % Residual
    res = muav - mul;
    % Newton increment
    a_inc = res/(asat*gaa);
    asat = asat + a_inc;
    %asat = min(asat,1);
    %asat = max(asat,0);
    wl = (ww - (1 - asat))/asat;
    
%     fprintf('Iter %.5g  residual %.5g  asat %.5g  wl %.5g \n',iter,res,asat,wl)
%     fprintf('muav %.5g  mul %.5g \n',muav,mul)
%     pause
    
end


end
