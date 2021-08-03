function [] = check_eos_sg(p0, T, eta0, q0, vareta, varq, constants)

format long

% Check derivatives of rho

% The quantity to check is the derivative of (rho' / rho^2) wrt p', eta', and q'

% This version takes account of subgrid variability in eta and q
% in calculating w-level quantities.  Cloud fraction and liquid water
% are calculated here.

% Useful constants
rr2 = sqrt(0.5);
rr2pi = sqrt(1/(2*pi));

% Assumed correlation between entropy and total water
% *** Put this in constants.param ***
rr = 0.0;

% Standard deviations
etastd = sqrt(vareta);
qstd   = sqrt(varq);


p = p0;
eta = eta0;
for q = q0 + [-1,1]*1e-4
    
    q
    
    % Quantities based on fluid means
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(p,T,q,constants.therm);
    for iter = 1:3
        res = eta + gt;
        T_inc = -res/gtt;
        T = T + T_inc;
        [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(p,T,q,constants.therm);
    end
    T
 
    drdp    = (gpt*gpt/gtt - gpp);
    drdeta  =  gpt/gtt;
    drdq    = (gpt*gtw/gtt - gpw)
    rhomean = 1/gp;
    qlmean = (a + q - 1)/a;
    
    % Determine density if all water were gas.
    % Need to iterate to find Tgas.
    Tgas = T;
    a = 1 - q;
    for iter = 1:3
        [g,gp,gt,ga,gpp,gpt,gpa,gtt,gta,gaa,gtwv] = gibbsav(p,Tgas,a,constants.therm);
        res = eta + gt;
        T_inc = -res/gtt;
        Tgas = Tgas + T_inc;
    end
    rhogas = 1/gp;
    drhogasdp = (gpt*gpt/gtt - gpp)*rhogas*rhogas;
    drhogasdeta = (gpt/gtt)*rhogas*rhogas;
    drhogasdq = (gpa - gpt*gta/gtt)*rhogas*rhogas;
    
    % Compute qsat and some derivatives at this eta
    [qsat,dqsatdeta,dqldq,dqldeta,ddrhodql,dqsatdp] = ...
           find_qsatl_new(p,eta,Tgas,constants.therm);
    
    % Difference between mean q and qsat
    Deltaq = q - qsat;
    
    % Alternative (better) estimate for dqldq and ddrhodql
    dqldq_alt = dqldq;
    ddrhodql_alt = ddrhodql; 
    if Deltaq > 0
        dqldq_alt = qlmean/Deltaq;
        ddrhodql_alt = (rhomean - rhogas)/qlmean;
    end
    
    % Standard deviation parameter
    sq2 = qstd*qstd ...
        - 2*rr*qstd*etastd*dqsatdeta ...
        + etastd*etastd*dqsatdeta*dqsatdeta;
    sq = sqrt(sq2);
    
    % Q1 parameter
    Q1 = Deltaq/sq;
    
    % Cloud fraction
    cldfrac = 0.5*(1 + erf(rr2*Q1));
    
    % Liquid water
    ql = dqldq_alt*(Deltaq*cldfrac + sq*rr2pi*exp(-0.5*(Deltaq/sq)^2));
    
    % Covariance of liquid and scalars
    Covarqleta = cldfrac*(dqldq_alt*qstd*etastd*rr + dqldeta*etastd*etastd);
    Covarqlq   = cldfrac*(dqldq_alt*qstd*qstd      + dqldeta*etastd*qstd*rr);
    
    % Liquid water variance
    Varql = sq2*(dqldq_alt^2)*cldfrac ...
               + dqldq_alt*Deltaq*ql ...
               - ql^2;
    
    % Density
    rhow = rhogas + ddrhodql_alt*ql
    
    % Coefficient needed for N^2
    ddrhodp = -ddrhodql_alt*dqldq_alt*cldfrac*dqsatdp;
    drho_parceldp = drhogasdp + ddrhodp;
    
    % Coefficients needed to compute correlation of scalars with buoyancy
    rho_deriv_eta = drhogasdeta + cldfrac*ddrhodql_alt*dqldeta;
    rho_deriv_q   = drhogasdq   + cldfrac*ddrhodql_alt*dqldq_alt;
    
    % Improved estimates of derivatives of density wrt p, eta, and q
    rrhosq = 1/rhow^2;
    drdpbar = drho_parceldp*rrhosq;
    drdeta  = rho_deriv_eta*rrhosq;
    drdq    = rho_deriv_q*rrhosq
    
    % Density variance
    Varrhow = (drhogasdeta*etastd)^2 ...
             + rr*drhogasdeta*drhogasdq*etastd*qstd ...
             + (drhogasdq*qstd)^2 ...
             + ddrhodql_alt*(dqldeta*Covarqleta + dqldq_alt*Covarqlq) ...
             + ddrhodql_alt^2*Varql;
        
                                    
end


format short
pause

