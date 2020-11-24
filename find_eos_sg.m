function eos = find_eos_sg(grid, state, constants)

% Determine densities and find residuals in equations of state

% Also compute terms needed in linearization. The quantity needed
% is the derivative of (rho' / rho^2) wrt p', eta', and q'

% This version takes account of subgrid variability in eta and q
% in calculating w-level quantities.  Cloud fraction and liquid water
% are calculated here.

% Useful constants
rr2 = sqrt(0.5);
rr2pi = sqrt(1/(2*pi));

% Assumed correlation between entropy and total water
% *** Put this in constants.param ***
rr = 0.0;

nz = grid.nz;
nzp = nz + 1;

% On p levels
for k = 1:nz
    
    % Fluid 1
    qbar   = grid.aboves(k)*state.fluid(1).q(k+1) ...
           + grid.belows(k)*state.fluid(1).q(k);
    etabar = grid.aboves(k)*state.fluid(1).eta(k+1) ...
           + grid.belows(k)*state.fluid(1).eta(k);
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(state.p(k),          ...
                                                   state.fluid(1).T(k), ...
                                                   qbar,                ...
                                                   constants.therm);
    eos.rho1(k) = 1.0/gp;
    eos.res_rho1(k) = 0.0;
    eos.sigma1(k) = state.fluid(1).m(k)*gp;
    eos.res_etap1(k) = etabar + gt;
    eos.drdp1(k)    = (gpt*gpt/gtt - gpp);
    eos.drdetap1(k) =  gpt/gtt;
    eos.drdqp1(k)   = (gpt*gtw/gtt - gpw);
    
    % Fluid 2
    qbar   = grid.aboves(k)*state.fluid(2).q(k+1) ...
           + grid.belows(k)*state.fluid(2).q(k);
    etabar = grid.aboves(k)*state.fluid(2).eta(k+1) ...
           + grid.belows(k)*state.fluid(2).eta(k);
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(state.p(k),          ...
                                                   state.fluid(2).T(k), ...
                                                   qbar,                ...
                                                   constants.therm);
    eos.rho2(k) = 1.0/gp;
    eos.res_rho2(k) = 0.0;
    eos.sigma2(k) = state.fluid(2).m(k)*gp;
    eos.res_etap2(k) = etabar + gt;
    eos.drdp2(k)    = (gpt*gpt/gtt - gpp);
    eos.drdetap2(k) =  gpt/gtt;
    eos.drdqp2(k)   = (gpt*gtw/gtt - gpw);

    % Residual in sigma equation
    eos.res_sigma(k) = 1 - eos.sigma1(k) - eos.sigma2(k);

end


% On w levels
for k = 1:nzp
    
    if k == 1
        pbar = grid.extrapb1*state.p(1) + grid.extrapb2*state.p(2);
    elseif k == nzp
        pbar = grid.extraptnz*state.p(nz) + grid.extraptnzm*state.p(nz-1);
    else
        pbar   = grid.abovew(k)*state.p(k) ...
               + grid.beloww(k)*state.p(k-1);
    end
    
    % Fluid 1
    % Quantities based on fluid means
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(pbar,                 ...
                                                   state.fluid(1).Tw(k), ...
                                                   state.fluid(1).q(k),  ...
                                                   constants.therm);
    eos.res_eta1(k) = state.fluid(1).eta(k) + gt;
    eos.drdpbar1(k) = (gpt*gpt/gtt - gpp);
    eos.drdeta1(k)  =  gpt/gtt;
    eos.drdq1(k)    = (gpt*gtw/gtt - gpw);
    rhomean = 1/gp;
    qlmean = (a + state.fluid(1).q(k) - 1)/a;
    
    % Determine density if all water were gas.
    % Need to iterate to find Tgas. [Better to include in main iteration?]
    Tgas = state.fluid(1).Tw(k);
    a = 1 - state.fluid(1).q(k);
    for iter = 1:3
        [g,gp,gt,ga,gpp,gpt,gpa,gtt,gta,gaa,gtwv] = gibbsav(pbar,Tgas,a,constants.therm);
        res = state.fluid(1).eta(k) + gt;
        T_inc = -res/gtt;
        Tgas = Tgas + T_inc;
    end
    rhogas = 1/gp;
    drhogasdp = (gpt*gpt/gtt - gpp)*rhogas*rhogas;
    drhogasdeta = (gpt/gtt)*rhogas*rhogas;
    drhogasdq = (gpt*gtw/gtt - gpw)*rhogas*rhogas;
    
    % Compute qsat and some derivatives at this eta
    [qsat,dqsatdeta,dqldq,dqldeta,ddrhodql,dqsatdp] = ...
           find_qsatl_new(pbar,state.fluid(1).eta(k),Tgas,constants.therm);
    
    % Difference between mean q and qsat
    Deltaq = state.fluid(1).q(k) - qsat;
    
    % Alternative (better) estimate for dqldq and ddrhodql
    dqldq_alt = dqldq;
    ddrhodql_alt = ddrhodql; 
    if Deltaq > 0
        dqldq_alt = qlmean/Deltaq;
        ddrhodql_alt = (rhomean - rhogas)/qlmean;
    end
    
     % Standard deviation parameter
    etastd = sqrt(state.fluid(1).vareta(k));
    qstd   = sqrt(state.fluid(1).varq(k));
    sq2 = qstd*qstd ...
        - 2*rr*qstd*etastd*dqsatdeta ...
        + etastd*etastd*dqsatdeta*dqsatdeta;
    sq = sqrt(sq2);
    
    % Q1 parameter
    Q1 = Deltaq/sq;
    
    % Cloud fraction
    eos.cldfrac1(k) = 0.5*(1 + erf(rr2*Q1));
    
    % Liquid water
    eos.ql1(k) = dqldq_alt*(Deltaq*eos.cldfrac1(k) + sq*rr2pi*exp(-0.5*(Deltaq/sq)^2));
    
    % Covariance of liquid and scalars
    eos.Covarqleta1(k) = eos.cldfrac1(k)*(dqldq_alt*qstd*etastd*rr + dqldeta*etastd*etastd);
    eos.Covarqlq1(k)   = eos.cldfrac1(k)*(dqldq_alt*qstd*qstd      + dqldeta*etastd*qstd*rr);
    
    % Liquid water variance
    eos.Varql1(k) = sq2*(dqldq_alt^2)*eos.cldfrac1(k) ...
                  + dqldq_alt*Deltaq*eos.ql1(k) ...
                  - eos.ql1(k)^2;
    
    % Density
    eos.rhow1(k) = rhogas + ddrhodql_alt*eos.ql1(k);
    
    % Coefficient needed for N^2
    ddrhodp = -ddrhodql_alt*dqldq_alt*eos.cldfrac1(k)*dqsatdp;
    drho_parceldp1(k) = drhogasdp + ddrhodp;
    
    % Coefficients needed to compute correlation of scalars with buoyancy
    eos.rho_deriv_eta1(k) = drhogasdeta + eos.cldfrac1(k)*ddrhodql_alt*dqldeta;
    eos.rho_deriv_q1(k)   = drhogasdq   + eos.cldfrac1(k)*ddrhodql_alt*dqldq_alt;
    
    % Density variance
    eos.Varrhow1(k) = (drhogasdeta*etastd)^2 ...
                    + rr*drhogasdeta*drhogasdq*etastd*qstd ...
                    + (drhogasdq*qstd)^2 ...
                    + ddrhodql_alt*(dqldeta*eos.Covarqleta1(k) + dqldq_alt*eos.Covarqlq1(k)) ...
                    + ddrhodql_alt^2*eos.Varql1(k);
    
    
    % Fluid 2
    % Quantities based on fluid means
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(pbar,                 ...
                                                   state.fluid(2).Tw(k), ...
                                                   state.fluid(2).q(k),  ...
                                                   constants.therm);
    eos.res_eta2(k) = state.fluid(2).eta(k) + gt;
    eos.drdpbar2(k) = (gpt*gpt/gtt - gpp);
    eos.drdeta2(k)  =  gpt/gtt;
    eos.drdq2(k)    = (gpt*gtw/gtt - gpw);
    rhomean = 1/gp;
    qlmean = (a + state.fluid(2).q(k) - 1)/a;
    
    % Determine density if all water were gas.
    % Need to iterate to find Tgas. [Better to include in main iteration?]
    Tgas = state.fluid(2).Tw(k);
    a = 1 - state.fluid(2).q(k);
    for iter = 1:3
        [g,gp,gt,ga,gpp,gpt,gpa,gtt,gta,gaa,gtwv] = gibbsav(pbar,Tgas,a,constants.therm);
        res = state.fluid(2).eta(k) + gt;
        T_inc = -res/gtt;
        Tgas = Tgas + T_inc;
    end
    rhogas = 1/gp;
    drhogasdp = (gpt*gpt/gtt - gpp)*rhogas*rhogas;
    drhogasdeta = (gpt/gtt)*rhogas*rhogas;
    drhogasdq = (gpt*gtw/gtt - gpw)*rhogas*rhogas;
    
    % Compute qsat and some derivatives at this eta
    [qsat,dqsatdeta,dqldq,dqldeta,ddrhodql,dqsatdp] = ...
           find_qsatl_new(pbar,state.fluid(2).eta(k),Tgas,constants.therm);
    
    % Difference between mean q and qsat
    Deltaq = state.fluid(2).q(k) - qsat;
    
    % Alternative (better) estimate for dqldq and ddrhodql
    dqldq_alt = dqldq;
    ddrhodql_alt = ddrhodql; 
    if Deltaq > 0
        dqldq_alt = qlmean/Deltaq;
        ddrhodql_alt = (rhomean - rhogas)/qlmean;
    end
    
     % Standard deviation parameter
    etastd = sqrt(state.fluid(2).vareta(k));
    qstd   = sqrt(state.fluid(2).varq(k));
    sq2 = qstd*qstd ...
        - 2*rr*qstd*etastd*dqsatdeta ...
        + etastd*etastd*dqsatdeta*dqsatdeta;
    sq = sqrt(sq2);
    
    % Q1 parameter
    Q1 = Deltaq/sq;
    
    % Cloud fraction
    eos.cldfrac2(k) = 0.5*(1 + erf(rr2*Q1));
    
    % Liquid water
    eos.ql2(k) = dqldq_alt*(Deltaq*eos.cldfrac2(k) + sq*rr2pi*exp(-0.5*(Deltaq/sq)^2));
        
    % Covariance of liquid and scalars
    eos.Covarqleta2(k) = eos.cldfrac2(k)*(dqldq_alt*qstd*etastd*rr + dqldeta*etastd*etastd);
    eos.Covarqlq2(k)   = eos.cldfrac2(k)*(dqldq_alt*qstd*qstd      + dqldeta*etastd*qstd*rr);
    
    % Liquid water variance
    eos.Varql2(k) = sq2*(dqldq_alt^2)*eos.cldfrac2(k) ...
                  + dqldq_alt*Deltaq*eos.ql2(k) ...
                  - eos.ql2(k)^2;
              
    % Density
    eos.rhow2(k) = rhogas + ddrhodql_alt*eos.ql2(k);

    % Coefficient needed for N^2
    ddrhodp = -ddrhodql_alt*dqldq_alt*eos.cldfrac2(k)*dqsatdp;
    drho_parceldp2(k) = drhogasdp + ddrhodp; 

    % Coefficients needed to compute correlation of scalars with buoyancy
    eos.rho_deriv_eta2(k) = drhogasdeta + eos.cldfrac2(k)*ddrhodql_alt*dqldeta;
    eos.rho_deriv_q2(k)   = drhogasdq   + eos.cldfrac2(k)*ddrhodql_alt*dqldq_alt;
    
    % Density variance
    eos.Varrhow2(k) = (drhogasdeta*etastd)^2 ...
                    + rr*drhogasdeta*drhogasdq*etastd*qstd ...
                    + (drhogasdq*qstd)^2 ...
                    + ddrhodql_alt*(dqldeta*eos.Covarqleta2(k) + dqldq_alt*eos.Covarqlq2(k)) ...
                    + ddrhodql_alt^2*eos.Varql2(k);
    
    eos.ccc1(k) = (drhogasdeta*etastd)^2 ...
                    + rr*drhogasdeta*drhogasdq*etastd*qstd ...
                    + (drhogasdq*qstd)^2;
    eos.ccc2(k) = ddrhodql_alt*(dqldeta*eos.Covarqleta2(k) + dqldq_alt*eos.Covarqlq2(k));
    eos.ccc3(k) = ddrhodql_alt^2*eos.Varql2(k);
    
    if eos.Varrhow2(k) < 0
        disp(['Negative rho variance  ',num2str(eos.Varrhow2(k)),'  k = ',num2str(k)])
        disp(['Gas = ',num2str(eos.ccc1(k)),...
              ' liquid = ',num2str(eos.ccc3(k)),...
              ' cross = ',num2str(eos.ccc2(k))])
        disp(['ql variance  ',num2str(eos.Varql2(k))])          
        disp([num2str(sq2*eos.cldfrac2(k)),'  ',...
              num2str(dqldq_alt*Deltaq*eos.ql2(k)),'  ',...
              num2str(- eos.ql2(k)^2)]);
        pause
    end
                
    eos.pbar(k) = pbar;
    
    % eos.theta1(k) = potential_temp(state.fluid(1).Tw(k),pbar,constants.phys,constants.therm);
    % eos.theta2(k) = potential_temp(state.fluid(2).Tw(k),pbar,constants.phys,constants.therm);
    eos.theta1(k) = eta2thetal( state.fluid(1).eta(k), ...
                                state.fluid(1).q(k), ...
                                state.fluid(1).Tw(k), ...
                                constants.therm, constants.phys.p00);
    eos.theta2(k) = eta2thetal( state.fluid(2).eta(k), ...
                                state.fluid(2).q(k), ...
                                state.fluid(2).Tw(k), ...
                                constants.therm, constants.phys.p00);
    eos.theta_rho1(k) = density_potential_temp(eos.rhow1(k),pbar,constants.phys,constants.therm);
    eos.theta_rho2(k) = density_potential_temp(eos.rhow2(k),pbar,constants.phys,constants.therm);

    eos.rho_ptl1(k) = potential_density(constants.phys.p00,state.fluid(1).eta(k),...
                                        state.fluid(1).q(k),state.fluid(1).Tw(k),constants.therm);
    eos.rho_ptl2(k) = potential_density(constants.phys.p00,state.fluid(2).eta(k),...
                                        state.fluid(2).q(k),state.fluid(2).Tw(k),constants.therm);
    
end


% Vertical pressure gradient
dpdz(2:nz)   = (state.p(2:nz) - state.p(1:nz-1))./grid.dzw(2:nz);
dpdz(1) = dpdz(2);
dpdz(nzp) = dpdz(nz);
%dpdzbar = grid.abovep.*dpdz(2:nzp) + grid.belowp.*dpdz(1:nz);

% Buoyancy frequency squared
% Note there are several ways to define this, depending on what is assumed
% about the parcel and the environment through which it moves. Here we
% consider a parcel that conserves its eta and q moving relative to a
% background in which eta1, eta2, q1, q2, and p are steady.

% Environmental density gradient
drho1dz = (eos.rhow1(2:nzp) - eos.rhow1(1:nz))./grid.dzp;
drho2dz = (eos.rhow2(2:nzp) - eos.rhow2(1:nz))./grid.dzp;
drhoenvdz = eos.sigma1.*drho1dz + eos.sigma2.*drho2dz;

% Parcel density gradient and N^2
temp = drho_parceldp1.*dpdz;
drho_parceldz = grid.abovep.*temp(2:nzp) + grid.belowp.*temp(1:nz);
eos.nsq1 = constants.phys.gravity*(drho_parceldz - drhoenvdz)./eos.rho1;
temp = drho_parceldp2.*dpdz;
drho_parceldz = grid.abovep.*temp(2:nzp) + grid.belowp.*temp(1:nz);
eos.nsq2 = constants.phys.gravity*(drho_parceldz - drhoenvdz)./eos.rho2;


%disp('** buoyancy correlation times 0.3 **')
%eos.rho_deriv_eta1(:) = 0.3*eos.rho_deriv_eta1(:);
%eos.rho_deriv_eta2(:) = 0.3*eos.rho_deriv_eta2(:);

