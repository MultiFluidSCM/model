function [ qsat, dqsatdeta ] = find_qsatl( p, eta, T, therm)
%FIND_QSATL Compute saturation mixing ratio.
%   Determine the mass mixing ratio of water needed to
%   make this air sample saturated at the given pressure p and
%   entropy eta. Also compute the derivative of qsat wrt eta.
%   A Newton method is used. T is a first guess for the
%   saturation temperature.
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

for iter = 1:3
    
    % Chemical potential (=Gibbs) of the liquid
    [mul,gp,glt,gpp,gpt,gtt] = gibbsl(p,T,therm);
    
    % Chemical potential of the vapour in the air-vapour mix
    [g,gp,gt,ga,gpp,gpt,gpa,gtt,gta,gaa,gtw] = gibbsav(p,T,asat,therm);
    muav = g - asat*ga;
    
    % Residuals
    res_mu = muav - mul;
    res_eta = gt + eta;
    
    % Coefficients in linearization
    R11 = gt - asat*gta - glt;
    R12 = -asat*gaa;
    R21 = gtt;
    R22 = gta;
    
    % Newton increment
    rdet = 1/(R11*R22 - R12*R21);
    T_inc = - rdet*(  R22*res_mu - R12*res_eta);
    a_inc = - rdet*(- R21*res_mu + R11*res_eta);
    T = T + T_inc;
    asat = asat + a_inc;
    
end

% On saturation curve ql = 0
qsat = 1 - asat;

% Derivative
Delta = gt - asat*gta - glt;
dqsatdeta = Delta/(Delta*gta + asat*gaa*gtt);


end
