function [ qsat, dqsatdeta, dqldq, dqldeta, ddrhodql, dqsatdp ] = find_qsatl_new( p, eta, T, therm)
%FIND_QSATL Compute saturation mixing ratio.
%   Determine the mass mixing ratio of water needed to
%   make this air sample saturated at the given pressure p and
%   entropy eta. Also compute the derivative of qsat wrt eta,
%   and the derivative of ql wrt q,
%   and the derivative of ql wrt eta,
%   and the derivative of rho - rhog wrt ql,
%   and the derivative of qsat wrt p.
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
    [mul,glp,glt,gpp,gpt,gtt] = gibbsl(p,T,therm);
    
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

% Derivatives
LambdaT = gt - asat*gta - glt;
Lambdap = gp - asat*gpa - glp;
den = LambdaT*LambdaT - asat*asat*gaa*gtt;
dqldq = -(asat*gaa*gtt + LambdaT*gta) / den;
dqldeta = LambdaT / den;
%dqsatdeta = LambdaT/(LambdaT*gta + asat*gaa*gtt);
dqsatdeta = - dqldeta/dqldq;
% dqldq = (asat*gaa*gtt + LambdaT*gta) / (asat*asat*gaa*gtt - LambdaT*LambdaT);
ddrhodql = -(LambdaT*gpt - Lambdap*gtt) / (gtt*gp*gp);
dqsatdp = (gpt*LambdaT - gtt*Lambdap)/(asat*gaa*gtt + LambdaT*gta);

end
