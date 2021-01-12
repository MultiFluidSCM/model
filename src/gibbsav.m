function [ g, gp, gt, ga, gpp, gpt, gpa, gtt, gta, gaa, gtwv ] = gibbsav( p, T, a, therm )
%GIBBSAV Gibbs function and derivatives for air - water vapour mix
%   Given pressure p, temperature T and mass fraction of dry air a,
%   compute the Gibbs thermodynamic potential for a mixture of dry
%   air and water vapour, along with its first and second derivatives.

%   This version assumes that dry air and water vapour are both perfect
%   gases and there are no virial terms representing interactions between
%   the two.

%   The latent heating term has been moved from vapour to liquid cf T17

% Unpack constants for clarity
epsilon = therm.epsilon;
T0      = therm.T0;
p0d     = therm.p0d;
p0v     = therm.p0v;
Cpd     = therm.Cpd;
Cpv     = therm.Cpv;
Rd      = therm.Rd;
Rv      = therm.Rv;


% Mass fraction of vapour in the mixture
b = 1 - a;

% Some useful intermediate quantities
den = 1 + (epsilon-1)*a;
LT = log(T/T0);
if (a == 0)
    L1 = 0;
    aL1a = 0;
    aL1aa = 0;
else
    L1 = log(epsilon*a*p/(den*p0d));
    aL1a =  L1 + 1/den;
    aL1aa = 1/(a*den*den);
end
if (b == 0)
    L2 = 0;
    bL2a = 0;
    bL2aa = 0;
else
    L2 = log(b*p/(den*p0v));
    bL2a = -L2 - epsilon/den;
    bL2aa = epsilon*epsilon/(b*den*den);
end

% Gibbs function
g = -(a*Cpd + b*Cpv)*LT*T ...
  + (a*Rd*L1 + b*Rv*L2)*T;

% First derivatives
gp = (a*Rd + b*Rv)*(T/p);
gt = -(a*Cpd + b*Cpv)*(1 + LT) +(a*Rd*L1 + b*Rv*L2);
ga = -(Cpd - Cpv)*LT*T ...
   + (Rd*aL1a + Rv*bL2a)*T;

% Second derivatives
gpp = -(a*Rd + b*Rv)*(T/(p*p));
gpt =  (a*Rd + b*Rv)/p;
gtt = -(a*Cpd + b*Cpv)/T;
gta = -(Cpd - Cpv)*(1 + LT) ...
    + (Rd*aL1a + Rv*bL2a);
gpa = (Rd - Rv)*(T/p);
gaa = (Rd*aL1aa + Rv*bL2aa)*T;

% -1 * specific entropy of water vapour
gtwv = - Cpv*(1 + LT) + Rv*L2;



end

