% Determine a correction to be added to residuals to mimic a converged
% solution

% The calculations here should mirror those in quasinewton



% ------

% Tendencies computed from new state
dt = time.dt;
t_new = time.t + dt;
old_diff.flag = 0;
[tend,relabel,eos,force,scales,surface_flux,budgets,work] = ...
    tendencies(grid,state_new,constants,t_new,dt,switches,old_diff);


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
tke1 = state_new.fluid(1).tke;
tke2 = state_new.fluid(2).tke;
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
    lhs1tke = state_new.fluid(2).tke.*m2    + state_new.fluid(1).tke.*m1    - adt*tend.fluid(1).mtke.tot;
    lhs2m   = - adt*tend.fluid(2).m.tot;
    lhs2eta = - adt*tend.fluid(2).meta.tot;
    lhs2q   = - adt*tend.fluid(2).mq.tot;
    lhs2w   = - adt*tend.fluid(2).mw.tot;
    lhs2u   = - adt*tend.fluid(2).mu.tot;
    lhs2v   = - adt*tend.fluid(2).mv.tot;
    lhs2tke = - adt*tend.fluid(2).mtke.tot;
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


% Save for later use
resfix1m   =  -res1m;
resfix1eta =  -res1eta;
resfix1q   =  -res1q;
resfix1w   =  -res1w;
resfix1u   =  -res1u;
resfix1v   =  -res1v;
resfix1tke =  -res1tke;
resfix2m   =  -res2m;
resfix2eta =  -res2eta;
resfix2q   =  -res2q;
resfix2w   =  -res2w;
resfix2u   =  -res2u;
resfix2v   =  -res2v;
resfix2tke =  -res2tke;
resfixeoseta1 = - eos.res_eta1;
resfixeoseta2 = - eos.res_eta2;
resfixeosetap1 = - eos.res_etap1;
resfixeosetap2 = - eos.res_etap2;
resfixeosrho1 = - eos.res_rho1;
resfixeosrho2 = - eos.res_rho2;
resfixeossigma = -eos.res_sigma;


% --------





