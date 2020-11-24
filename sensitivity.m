% For diagnostics.
% Examine sensitivity of residuals to changes in unknowns


% Levels for displaying diagnostics
krange_sens = 51:53;

% ------

% Save state before perturbing
state_save = state_new;

% ------

% Define perturbations
pp_p    = zeros(1,nz);
pp_m1   = zeros(1,nz);
pp_w1   = zeros(1,nzp);
pp_eta1 = zeros(1,nzp);
pp_q1   = zeros(1,nzp);
pp_T1   = zeros(1,nz);
pp_Tw1  = zeros(1,nzp);
pp_m2   = zeros(1,nz);
pp_w2   = zeros(1,nzp);
pp_eta2 = zeros(1,nzp);
pp_q2   = zeros(1,nzp);
pp_T2   = zeros(1,nz);
pp_Tw2  = zeros(1,nzp);
pp_tke1 = zeros(1,nz);
pp_tke2 = zeros(1,nz);
pp_vareta1 = zeros(1,nzp);
pp_vareta2 = zeros(1,nzp);
pp_varq1 = zeros(1,nzp);
pp_varq2 = zeros(1,nzp);

pp_p(56) = 0.0;
pp_m1(56) =  0.0;
pp_m2(56) =  0.0;
pp_w1(56) =  0.0;
pp_w2(56) =  0.0;
pp_eta1(56) =  0.0;
pp_eta2(56) =  0.0;
pp_q1(56) =  0.0000;
pp_q2(56) =  0.0000;
pp_tke1(56) = 0.0;
pp_tke2(28) = 0.0;
pp_vareta1(52) = 0.0;
pp_vareta2(52) = 0.1;

% Pack into one array
xx(1:9:9*nz+1) = pp_w1;
xx(2:9:9*nz+2) = pp_w2;
xx(3:9:9*nz+3) = pp_eta1;
xx(4:9:9*nz+4) = pp_eta2;
xx(5:9:9*nz+5) = pp_q1;
xx(6:9:9*nz+6) = pp_q2;
xx(7:9:9*nz-2) = pp_m1;
xx(8:9:9*nz-1) = pp_m2;
xx(9:9:9*nz)   = pp_p;

% ------

% Tendencies 

% Tendencies computed from new state
dt = time.dt;
t_new = time.t + dt;
[tend,relabel,eos,force,scales,surface_flux,budgets,work] = tendencies(grid,state_new,constants,t_new,dt,switches);

% Unpack some fields for clarity of code
m1 = state_new.fluid(1).m;
m2 = state_new.fluid(2).m;
w1 = state_new.fluid(1).w;
w2 = state_new.fluid(2).w;
eta1 = state_new.fluid(1).eta;
eta2 = state_new.fluid(2).eta;
q1 = state_new.fluid(1).q;
q2 = state_new.fluid(2).q;
u1 = state_new.fluid(1).u;
u2 = state_new.fluid(2).u;
v1 = state_new.fluid(1).v;
v2 = state_new.fluid(2).v;
m1bar = work.m1bar;
m2bar = work.m2bar;
dpdz = work.dpdz;
M12 = relabel.M12;
M21 = relabel.M21;
M12bar = relabel.M12bar;
M21bar = relabel.M21bar;
what12_w1     = relabel.what12   - work.w1ubar;
what21_w1     = relabel.what21   - work.w1ubar;
what12_w2     = relabel.what12   - work.w2ubar;
what21_w2     = relabel.what21   - work.w2ubar;
etahat12_eta1 = relabel.etahat12 - work.eta1ubar;
etahat21_eta1 = relabel.etahat21 - work.eta1ubar;
etahat12_eta2 = relabel.etahat12 - work.eta2ubar;
etahat21_eta2 = relabel.etahat21 - work.eta2ubar;
qhat12_q1     = relabel.qhat12   - work.q1ubar;
qhat21_q1     = relabel.qhat21   - work.q1ubar;
qhat12_q2     = relabel.qhat12   - work.q2ubar;
qhat21_q2     = relabel.qhat21   - work.q2ubar;

% Left hand sides of all equations
adt = time.alpha*time.dt;
if switches.c
    lhs1m   = m2                            + m1                            - adt*tend.fluid(1).m.tot;
    lhs1eta = state_new.fluid(2).eta.*m2bar + state_new.fluid(1).eta.*m1bar - adt*tend.fluid(1).meta.tot;
    lhs1q   = state_new.fluid(2).q  .*m2bar + state_new.fluid(1).q  .*m1bar - adt*tend.fluid(1).mq.tot;
    lhs1w   = state_new.fluid(2).w  .*m2bar + state_new.fluid(1).w  .*m1bar - adt*tend.fluid(1).mw.tot;
    lhs1u   = state_new.fluid(2).u  .*m2    + state_new.fluid(1).u  .*m1    - adt*tend.fluid(1).mu.tot;
    lhs1v   = state_new.fluid(2).v  .*m2    + state_new.fluid(1).v  .*m1    - adt*tend.fluid(1).mv.tot;
    lhs2m   = - adt*tend.fluid(2).m.tot;
    lhs2eta = - adt*tend.fluid(2).meta.tot;
    lhs2q   = - adt*tend.fluid(2).mq.tot;
    lhs2w   = - adt*tend.fluid(2).mw.tot;
    lhs2u   = - adt*tend.fluid(2).mu.tot;
    lhs2v   = - adt*tend.fluid(2).mv.tot;
else
    lhs1m   = m1                            - adt*tend.fluid(1).m.tot;
    lhs1eta = state_new.fluid(1).eta.*m1bar - adt*tend.fluid(1).meta.tot;
    lhs1q   = state_new.fluid(1).q  .*m1bar - adt*tend.fluid(1).mq.tot;
    lhs1w   = state_new.fluid(1).w  .*m1bar - adt*tend.fluid(1).mw.tot;
    lhs1u   = state_new.fluid(1).u  .*m1    - adt*tend.fluid(1).mu.tot;
    lhs1v   = state_new.fluid(1).v  .*m1    - adt*tend.fluid(1).mv.tot;
    lhs1tke = state_new.fluid(1).tke.*m1    - adt*tend.fluid(1).mtke.tot;
    lhs2m   = m2                            - adt*tend.fluid(2).m.tot;
    lhs2eta = state_new.fluid(2).eta.*m2bar - adt*tend.fluid(2).meta.tot;
    lhs2q   = state_new.fluid(2).q  .*m2bar - adt*tend.fluid(2).mq.tot;
    lhs2w   = state_new.fluid(2).w  .*m2bar - adt*tend.fluid(2).mw.tot;
    lhs2u   = state_new.fluid(2).u  .*m2    - adt*tend.fluid(2).mu.tot;
    lhs2v   = state_new.fluid(2).v  .*m2    - adt*tend.fluid(2).mv.tot;
    lhs2tke = state_new.fluid(2).tke.*m2    - adt*tend.fluid(2).mtke.tot;
end

% Residuals in all equations (note opposite sign convention
% to Gibbs paper)
res1m   =  rhs1m   - lhs1m;
res1eta =  rhs1eta - lhs1eta;
res1q   =  rhs1q   - lhs1q;
res1w   =  rhs1w   - lhs1w;
res1u   =  rhs1u   - lhs1u;
res1v   =  rhs1v   - lhs1v;
res1tke =  rhs1tke - lhs1tke;
res2m   =  rhs2m   - lhs2m;
res2eta =  rhs2eta - lhs2eta;
res2q   =  rhs2q   - lhs2q;
res2w   =  rhs2w   - lhs2w;
res2u   =  rhs2u   - lhs2u;
res2v   =  rhs2v   - lhs2v;
res2tke =  rhs2tke - lhs2tke;

% Interpolate mass residuals to w levels
% using `reversed' weighting for conservation
res1mbar = weight_to_w(grid,res1m);
res2mbar = weight_to_w(grid,res2m);

% Residuals in w equation after substituting for buoyancy
res1w = res1w + adt*m1bar.*dpdz.*eos.drdeta1.*eos.res_eta1;
res2w = res2w + adt*m2bar.*dpdz.*eos.drdeta2.*eos.res_eta2;

% Implied residuals in quasi-advective form eta, q, w, u and v equations
res1eta = res1eta - work.eta1ubar.*res1mbar;
res2eta = res2eta - work.eta2ubar.*res2mbar;
res1q   = res1q   - work.q1ubar  .*res1mbar;
res2q   = res2q   - work.q2ubar  .*res2mbar;
res1w   = res1w   - work.w1ubar  .*res1mbar;
res2w   = res2w   - work.w2ubar  .*res2mbar;
res1u   = res1u   - work.u1ubar  .*res1m;
res2u   = res2u   - work.u2ubar  .*res2m;
res1v   = res1v   - work.v1ubar  .*res1m;
res2v   = res2v   - work.v2ubar  .*res2m;

% Residual in sigma equation after substituting for temperature
res_s = eos.res_sigma + m1.*(eos.res_rho1 + eos.drdetap1.*eos.res_etap1) ...
                      + m2.*(eos.res_rho2 + eos.drdetap2.*eos.res_etap2);

disp(' ')
disp('Residuals')
% disp(['res_w1    ' num2str(res1w(krange_sens))])
% disp(['res_w2    ' num2str(res2w(krange_sens))])
% disp(['res_m1    ' num2str(res1m(krange_sens))])
% disp(['res_m2    ' num2str(res2m(krange_sens))])
%disp(['res_eta1  ' num2str(res1eta(krange_sens))])
%disp(['res_eta2  ' num2str(res2eta(krange_sens))])
% disp(['res_q1    ' num2str(res1q(krange_sens))])
% disp(['res_q2    ' num2str(res2q(krange_sens))])
%disp(['res_u1    ' num2str(res1u(krange_sens))])
%disp(['res_u2    ' num2str(res2u(krange_sens))])
%disp(['res_v1    ' num2str(res1v(krange_sens))])
%disp(['res_v2    ' num2str(res2v(krange_sens))])
%disp(['res_sigma ' num2str(res_s(krange_sens))])
%disp(['res_etap1 ' num2str(eos.res_etap1(krange_sens))])
%disp(['res_eta1  ' num2str(eos.res_eta1(krange_sens))])
%disp(['res_etap2 ' num2str(eos.res_etap2(krange_sens))])
%disp(['res_eta2  ' num2str(eos.res_eta2(krange_sens))])
%disp(['res_tke1  ' num2str(res1tke(krange_sens))])
%disp(['res_tke2  ' num2str(res2tke(krange_sens))])
disp(['tend_mvareta1  ' num2str(tend.fluid(1).mvareta.tot(krange_sens))])
disp(['tend_mvareta2  ' num2str(tend.fluid(2).mvareta.tot(krange_sens))])


% Save for later
oldtend = tend;
oldres_s   = res_s;
oldres1m   = res1m;
oldres2m   = res2m;
oldres1eta = res1eta;
oldres2eta = res2eta;
oldres1q   = res1q;
oldres2q   = res2q;
oldres1w   = res1w;
oldres2w   = res2w;
oldres1tke = res1tke;
oldres2tke = res2tke;

% ------

% Now perturb all variables
state_new.p            =     state_new.p            + pp_p;
state_new.fluid(1).m   =     state_new.fluid(1).m   + pp_m1;
state_new.fluid(1).w   =     state_new.fluid(1).w   + pp_w1;
state_new.fluid(1).eta =     state_new.fluid(1).eta + pp_eta1;
state_new.fluid(1).q   = max(state_new.fluid(1).q   + pp_q1, 0);
state_new.fluid(1).T   =     state_new.fluid(1).T   + pp_T1;
state_new.fluid(1).Tw  =     state_new.fluid(1).Tw  + pp_Tw1;
state_new.fluid(1).tke =     state_new.fluid(1).tke + pp_tke1;
state_new.fluid(1).vareta =  state_new.fluid(1).vareta + pp_vareta1;
state_new.fluid(1).varq   =  state_new.fluid(1).varq   + pp_varq1;
state_new.fluid(2).m   =     state_new.fluid(2).m   + pp_m2;
state_new.fluid(2).w   =     state_new.fluid(2).w   + pp_w2;
state_new.fluid(2).eta =     state_new.fluid(2).eta + pp_eta2;
state_new.fluid(2).q   = max(state_new.fluid(2).q   + pp_q2, 0);
state_new.fluid(2).T   =     state_new.fluid(2).T   + pp_T2;
state_new.fluid(2).Tw  =     state_new.fluid(2).Tw  + pp_Tw2;
state_new.fluid(2).tke =     state_new.fluid(2).tke + pp_tke2;
state_new.fluid(2).vareta =  state_new.fluid(2).vareta + pp_vareta2;
state_new.fluid(2).varq   =  state_new.fluid(2).varq   + pp_varq2;

% ------

% and repeat the calculation of residuals

% Tendencies computed from new state
dt = time.dt;
t_new = time.t + dt;
[tend,relabel,eos,force,scales,surface_flux,budgets,work] = tendencies(grid,state_new,constants,t_new,dt,switches);

% Unpack some fields for clarity of code
m1 = state_new.fluid(1).m;
m2 = state_new.fluid(2).m;
w1 = state_new.fluid(1).w;
w2 = state_new.fluid(2).w;
eta1 = state_new.fluid(1).eta;
eta2 = state_new.fluid(2).eta;
q1 = state_new.fluid(1).q;
q2 = state_new.fluid(2).q;
u1 = state_new.fluid(1).u;
u2 = state_new.fluid(2).u;
v1 = state_new.fluid(1).v;
v2 = state_new.fluid(2).v;
m1bar = work.m1bar;
m2bar = work.m2bar;
dpdz = work.dpdz;
M12 = relabel.M12;
M21 = relabel.M21;
M12bar = relabel.M12bar;
M21bar = relabel.M21bar;
what12_w1     = relabel.what12   - work.w1ubar;
what21_w1     = relabel.what21   - work.w1ubar;
what12_w2     = relabel.what12   - work.w2ubar;
what21_w2     = relabel.what21   - work.w2ubar;
etahat12_eta1 = relabel.etahat12 - work.eta1ubar;
etahat21_eta1 = relabel.etahat21 - work.eta1ubar;
etahat12_eta2 = relabel.etahat12 - work.eta2ubar;
etahat21_eta2 = relabel.etahat21 - work.eta2ubar;
qhat12_q1     = relabel.qhat12   - work.q1ubar;
qhat21_q1     = relabel.qhat21   - work.q1ubar;
qhat12_q2     = relabel.qhat12   - work.q2ubar;
qhat21_q2     = relabel.qhat21   - work.q2ubar;

% Left hand sides of all equations
adt = time.alpha*time.dt;
if switches.c
    lhs1m   = m2                            + m1                            - adt*tend.fluid(1).m.tot;
    lhs1eta = state_new.fluid(2).eta.*m2bar + state_new.fluid(1).eta.*m1bar - adt*tend.fluid(1).meta.tot;
    lhs1q   = state_new.fluid(2).q  .*m2bar + state_new.fluid(1).q  .*m1bar - adt*tend.fluid(1).mq.tot;
    lhs1w   = state_new.fluid(2).w  .*m2bar + state_new.fluid(1).w  .*m1bar - adt*tend.fluid(1).mw.tot;
    lhs1u   = state_new.fluid(2).u  .*m2    + state_new.fluid(1).u  .*m1    - adt*tend.fluid(1).mu.tot;
    lhs1v   = state_new.fluid(2).v  .*m2    + state_new.fluid(1).v  .*m1    - adt*tend.fluid(1).mv.tot;
    lhs2m   = - adt*tend.fluid(2).m.tot;
    lhs2eta = - adt*tend.fluid(2).meta.tot;
    lhs2q   = - adt*tend.fluid(2).mq.tot;
    lhs2w   = - adt*tend.fluid(2).mw.tot;
    lhs2u   = - adt*tend.fluid(2).mu.tot;
    lhs2v   = - adt*tend.fluid(2).mv.tot;
else
    lhs1m   = m1                            - adt*tend.fluid(1).m.tot;
    lhs1eta = state_new.fluid(1).eta.*m1bar - adt*tend.fluid(1).meta.tot;
    lhs1q   = state_new.fluid(1).q  .*m1bar - adt*tend.fluid(1).mq.tot;
    lhs1w   = state_new.fluid(1).w  .*m1bar - adt*tend.fluid(1).mw.tot;
    lhs1u   = state_new.fluid(1).u  .*m1    - adt*tend.fluid(1).mu.tot;
    lhs1v   = state_new.fluid(1).v  .*m1    - adt*tend.fluid(1).mv.tot;
    lhs1tke = state_new.fluid(1).tke.*m1    - adt*tend.fluid(1).mtke.tot;
    lhs2m   = m2                            - adt*tend.fluid(2).m.tot;
    lhs2eta = state_new.fluid(2).eta.*m2bar - adt*tend.fluid(2).meta.tot;
    lhs2q   = state_new.fluid(2).q  .*m2bar - adt*tend.fluid(2).mq.tot;
    lhs2w   = state_new.fluid(2).w  .*m2bar - adt*tend.fluid(2).mw.tot;
    lhs2u   = state_new.fluid(2).u  .*m2    - adt*tend.fluid(2).mu.tot;
    lhs2v   = state_new.fluid(2).v  .*m2    - adt*tend.fluid(2).mv.tot;
    lhs2tke = state_new.fluid(2).tke.*m2    - adt*tend.fluid(2).mtke.tot;
end

% Residuals in all equations (note opposite sign convention
% to Gibbs paper)
res1m   =  rhs1m   - lhs1m;
res1eta =  rhs1eta - lhs1eta;
res1q   =  rhs1q   - lhs1q;
res1w   =  rhs1w   - lhs1w;
res1u   =  rhs1u   - lhs1u;
res1v   =  rhs1v   - lhs1v;
res1tke =  rhs1tke - lhs1tke;
res2m   =  rhs2m   - lhs2m;
res2eta =  rhs2eta - lhs2eta;
res2q   =  rhs2q   - lhs2q;
res2w   =  rhs2w   - lhs2w;
res2u   =  rhs2u   - lhs2u;
res2v   =  rhs2v   - lhs2v;
res2tke =  rhs2tke - lhs2tke;

% Interpolate mass residuals to w levels
% using `reversed' weighting for conservation
res1mbar = weight_to_w(grid,res1m);
res2mbar = weight_to_w(grid,res2m);

% Residuals in w equation after substituting for buoyancy
res1w = res1w + adt*m1bar.*dpdz.*eos.drdeta1.*eos.res_eta1;
res2w = res2w + adt*m2bar.*dpdz.*eos.drdeta2.*eos.res_eta2;

% Implied residuals in quasi-advective form eta, q, w, u and v equations
res1eta = res1eta - work.eta1ubar.*res1mbar;
res2eta = res2eta - work.eta2ubar.*res2mbar;
res1q   = res1q   - work.q1ubar  .*res1mbar;
res2q   = res2q   - work.q2ubar  .*res2mbar;
res1w   = res1w   - work.w1ubar  .*res1mbar;
res2w   = res2w   - work.w2ubar  .*res2mbar;
res1u   = res1u   - work.u1ubar  .*res1m;
res2u   = res2u   - work.u2ubar  .*res2m;
res1v   = res1v   - work.v1ubar  .*res1m;
res2v   = res2v   - work.v2ubar  .*res2m;

% Residual in sigma equation after substituting for temperature
res_s = eos.res_sigma + m1.*(eos.res_rho1 + eos.drdetap1.*eos.res_etap1) ...
                      + m2.*(eos.res_rho2 + eos.drdetap2.*eos.res_etap2);

disp(' ')
disp('Residuals')
% disp(['res_w1    ' num2str(res1w(krange_sens))])
% disp(['res_w2    ' num2str(res2w(krange_sens))])
% disp(['res_m1    ' num2str(res1m(krange_sens))])
% disp(['res_m2    ' num2str(res2m(krange_sens))])
%disp(['res_eta1  ' num2str(res1eta(krange_sens))])
%disp(['res_eta2  ' num2str(res2eta(krange_sens))])
%disp(['res_q1    ' num2str(res1q(krange_sens))])
%disp(['res_q2    ' num2str(res2q(krange_sens))])
%disp(['res_u1    ' num2str(res1u(krange_sens))])
%disp(['res_u2    ' num2str(res2u(krange_sens))])
%disp(['res_v1    ' num2str(res1v(krange_sens))])
%disp(['res_v2    ' num2str(res2v(krange_sens))])
%disp(['res_sigma ' num2str(res_s(krange_sens))])
%disp(['res_etap1 ' num2str(eos.res_etap1(krange_sens))])
%disp(['res_eta1  ' num2str(eos.res_eta1(krange_sens))])
%disp(['res_etap2 ' num2str(eos.res_etap2(krange_sens))])
%disp(['res_eta2  ' num2str(eos.res_eta2(krange_sens))])
%disp(['res_tke1  ' num2str(res1tke(krange_sens))])
%disp(['res_tke2  ' num2str(res2tke(krange_sens))])
disp(['tend_mvareta1  ' num2str(tend.fluid(1).mvareta.tot(krange_sens))])
disp(['tend_mvareta2  ' num2str(tend.fluid(2).mvareta.tot(krange_sens))])

disp(' ')
disp('Change in residuals')
% disp(['res_w1    ' num2str(res1w(krange_sens)-oldres1w(krange_sens))])
% disp(['res_w2    ' num2str(res2w(krange_sens)-oldres2w(krange_sens))])
% disp(['res_m1    ' num2str(res1m(krange_sens)-oldres1m(krange_sens))])
% disp(['res_m2    ' num2str(res2m(krange_sens)-oldres2m(krange_sens))])
%disp(['res_eta1  ' num2str(res1eta(krange_sens)-oldres1eta(krange_sens))])
%disp(['res_eta2  ' num2str(res2eta(krange_sens)-oldres2eta(krange_sens))])
% disp(['res_q1    ' num2str(res1q(krange_sens)-oldres1q(krange_sens))])
% disp(['res_q2    ' num2str(res2q(krange_sens)-oldres2q(krange_sens))])
%disp(['res_u1    ' num2str(res1u(krange_sens))])
%disp(['res_u2    ' num2str(res2u(krange_sens))])
%disp(['res_v1    ' num2str(res1v(krange_sens))])
%disp(['res_v2    ' num2str(res2v(krange_sens))])
%disp(['res_sigma ' num2str(res_s(krange_sens)-oldres_s(krange_sens))])
%disp(['res_etap1 ' num2str(eos.res_etap1(krange_sens))])
%disp(['res_eta1  ' num2str(eos.res_eta1(krange_sens))])
%disp(['res_etap2 ' num2str(eos.res_etap2(krange_sens))])
%disp(['res_eta2  ' num2str(eos.res_eta2(krange_sens))])
%disp(['res_tke1    ' num2str(res1tke(krange_sens)-oldres1tke(krange_sens))])
%disp(['res_tke2    ' num2str(res2tke(krange_sens)-oldres2tke(krange_sens))])
disp(['tend_mvareta1  ' num2str(tend.fluid(1).mvareta.tot(krange_sens)-oldtend.fluid(1).mvareta.tot(krange_sens))])
disp(['tend_mvareta2  ' num2str(tend.fluid(2).mvareta.tot(krange_sens)-oldtend.fluid(2).mvareta.tot(krange_sens))])

% --------

% Build linear system for w-eta-q-m-p system
build_linear_system

% ---------

% For comparison, compute the predicted change in residuals
rr = - Ndiagmult(cc,xx);
rr_w1   = rr(1:9:9*nz+1);
rr_w2   = rr(2:9:9*nz+2);
rr_eta1 = rr(3:9:9*nz+3);
rr_eta2 = rr(4:9:9*nz+4);
rr_q1   = rr(5:9:9*nz+5);
rr_q2   = rr(6:9:9*nz+6);
rr_m1   = rr(7:9:9*nz-2);
rr_m2   = rr(8:9:9*nz-1);
rr_s    = rr(9:9:9*nz);

disp(' ')
disp('Predicted change in residuals')
% disp(['res_w1    ' num2str(rr_w1(krange_sens))])
% disp(['res_w2    ' num2str(rr_w2(krange_sens))])
% disp(['res_m1    ' num2str(rr_m1(krange_sens))])
% disp(['res_m2    ' num2str(rr_m2(krange_sens))])
% disp(['res_eta1  ' num2str(rr_eta1(krange_sens))])
% disp(['res_eta2  ' num2str(rr_eta2(krange_sens))])
% disp(['res_q1    ' num2str(rr_q1(krange_sens))])
% disp(['res_q2    ' num2str(rr_q2(krange_sens))])
% disp(['res_sigma ' num2str(rr_s(krange_sens))])


% Build linear system for tke1 system
ix = 1:nz;
dd = zeros(3,nz);
% Tendency term
dd(2,ix) = dd(2,ix) + m1;
% Transport terms    
dd(1,ix) = dd(1,ix) - adt*work.dFtke1dtkeb(1:nz)./dzp;
dd(2,ix) = dd(2,ix) - adt*(work.dFtke1dtkea(1:nz) - work.dFtke1dtkeb(2:nzp))./dzp;
dd(3,ix) = dd(3,ix) + adt*work.dFtke1dtkea(2:nzp)./dzp;
% Diffusion terms
dd(1,ix) = dd(1,ix) - adt*work.dDtke1dtkeb(1:nz)./dzp;
dd(2,ix) = dd(2,ix) - adt*(work.dDtke1dtkea(1:nz) - work.dDtke1dtkeb(2:nzp))./dzp;
dd(3,ix) = dd(3,ix) + adt*work.dDtke1dtkea(2:nzp)./dzp;
% Dissipation term
%dd(2,ix) = dd(2,ix) + adt*m1.*(1.5./scales.L_turb1 - tke1.*work.dLdtke1./(scales.L_turb1.^2)).*sqrt(tke1);
dd(2,ix) = dd(2,ix) + adt*m1.*(1.5./scales.L_turb1).*sqrt(tke1);
% Relabelling terms - only keep dependence on tke1
dd(2,ix) = dd(2,ix) + adt*M21;

% Predicted residual in tke1 equation
rr_tke1 = - Ndiagmult(dd,pp_tke1);
%disp(['res_tke1 ' num2str(rr_tke1(krange_sens))])


% Build linear system for tke2 system
ix = 1:nz;
dd = zeros(3,nz);
% Tendency term
dd(2,ix) = dd(2,ix) + m2;
% Transport terms    
dd(1,ix) = dd(1,ix) - adt*work.dFtke2dtkeb(1:nz)./dzp;
dd(2,ix) = dd(2,ix) - adt*(work.dFtke2dtkea(1:nz) - work.dFtke2dtkeb(2:nzp))./dzp;
dd(3,ix) = dd(3,ix) + adt*work.dFtke2dtkea(2:nzp)./dzp;
% Diffusion terms
dd(1,ix) = dd(1,ix) - adt*work.dDtke2dtkeb(1:nz)./dzp;
dd(2,ix) = dd(2,ix) - adt*(work.dDtke2dtkea(1:nz) - work.dDtke2dtkeb(2:nzp))./dzp;
dd(3,ix) = dd(3,ix) + adt*work.dDtke2dtkea(2:nzp)./dzp;
% Dissipation term
%dd(2,ix) = dd(2,ix) + adt*m2.*(1.5./scales.L_turb2 - tke2.*work.dLdtke2./(scales.L_turb2.^2)).*sqrt(tke2);
dd(2,ix) = dd(2,ix) + adt*m2.*(1.5./scales.L_turb2).*sqrt(tke2);
% Relabelling terms - only keep dependence on tke2
dd(2,ix) = dd(2,ix) + adt*M12;

% Predicted residual in tke1 equation
rr_tke2 = - Ndiagmult(dd,pp_tke2);
%disp(['res_tke2 ' num2str(rr_tke2(krange_sens))])



% Linear systems for variance equations

% Turbulence time scales at w-levels
T_turb1_bar = weight_to_w(grid,scales.T_turb1);
T_turb2_bar = weight_to_w(grid,scales.T_turb2);

% and for dissipation of flux in buoyancy correlation terms
t_scale1 = 6*scales.L_turb1./sqrt(tke1);
t_scale2 = 6*scales.L_turb2./sqrt(tke2);

% Factors needed to allow for buoyancy correlation terms in
% linearization
deta1dz = (eta1(2:nzp) - eta1(1:nz))./grid.dzp;
deta2dz = (eta2(2:nzp) - eta2(1:nz))./grid.dzp;
temp = constants.phys.gravity*t_scale1.*deta1dz.*eos.sigma1;
factor1(2:nzp) = belowr(2:nzp).*abovep.*temp;
factor1(1) = 0;
factor1(1:nz) = factor1(1:nz) + abover(1:nz).*belowp.*temp;
factor1 = -factor1.*eos.rho_deriv_eta1;
temp = constants.phys.gravity*t_scale2.*deta2dz.*eos.sigma2;
factor2(2:nzp) = belowr(2:nzp).*abovep.*temp;
factor2(1) = 0;
factor2(1:nz) = factor2(1:nz) + abover(1:nz).*belowp.*temp;
factor2 = -factor2.*eos.rho_deriv_eta2;

for k = 1:nzp

    % Set up 2x2 linear system for eta variance
    detfac = M12bar(k)*relabel.f_sort_chi_hat(k)/max(0.001,sqrt(state_new.fluid(2).vareta(k)));
    A11 =  M12bar(k) + m1bar(k)/T_turb1_bar(k); % + factor1(k);
    A12 = -M12bar(k) - (relabel.etahat12(k) - eta1(k))*detfac;
    A21 = -M21bar(k);
    A22 =  M21bar(k) + m2bar(k)/T_turb2_bar(k) + (relabel.etahat12(k) - eta2(k))*detfac; % + factor2(k);
    
    % And solve the linear system
    rr_vareta1(k) = -(A11*pp_vareta1(k) + A12*pp_vareta2(k));
    rr_vareta2(k) = -(A21*pp_vareta1(k) + A22*pp_vareta2(k));
                     
    % Set up 2x2 linear system for q variance
    detfac = M12bar(k)*relabel.f_sort_chi_hat(k)/max(1e-6,sqrt(state_new.fluid(2).varq(k)));
    A11 =  M12bar(k) + m1bar(k)/T_turb1_bar(k);
    A12 = -M12bar(k) - (relabel.qhat12(k) - q1(k))*detfac;
    A21 = -M21bar(k);
    A22 =  M21bar(k) + m2bar(k)/T_turb2_bar(k) + (relabel.qhat12(k) - q2(k))*detfac;

    % And solve the linear system
    rr_varq1(k) = -(A11*pp_varq1(k) + A12*pp_varq2(k));
    rr_varq2(k) = -(A21*pp_varq1(k) + A22*pp_varq2(k));

end

disp(['tend_mvareta1  ' num2str(rr_vareta1(krange_sens))])
disp(['tend_mvareta2  ' num2str(rr_vareta2(krange_sens))])


% Alternatively
build_etavar_system

pp_var(1:2:2*nzp-1) = pp_vareta1;
pp_var(2:2:2*nzp  ) = pp_vareta2;
rr = -Ndiagmult(cc,pp_var);
rr_vareta1 = rr(1:2:2*nzp-1);
rr_vareta2 = rr(2:2:2*nzp  );
disp('or...')
disp(['tend_mvareta1  ' num2str(rr_vareta1(krange_sens))])
disp(['tend_mvareta2  ' num2str(rr_vareta2(krange_sens))])


pause

% --------

% Restore state_new from saved state
state_new = state_save;

% --------


