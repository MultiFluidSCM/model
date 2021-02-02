% Advance the integration by one step
% Assume that the first guess for the new
% time level solution has already been set to
% the old time level solution

disp(' ')
disp(['Step ',num2str(istep),' of ',num2str(time.nstop),'(',num2str(100*istep/time.nstop),'%)'])

% if (istep == 5)
%     % To check equal tendencies result from equal fluid properties
%     % call equate_fluid_properties here and
%     % call check_equal_tendencies from tendencies
%     disp('*** Equating fluid properties ***')
%     equate_fluid_properties
% end

% ------

% Unpack some fields needed below
m1   = state_new.fluid(1).m;
m2   = state_new.fluid(2).m;

% ------

% Tendencies computed from current state
old_diff.flag = 0;
[tend,relabel,eos,force,scales,surface_flux,budgets,work] = ...
         tendencies(grid,state_old,constants,time.t,time.dt,switches,old_diff);
m1bar = work.m1bar;
m2bar = work.m2bar;

% Save tendencies and diffusivities for diagnostics/testing
% disp('*** Saving old tendencies ***')
% tend_old = tend;
% disp('*** Saving old diffusivities ***')
% old_diff.kdifft1 = work.kdifft1;
% old_diff.kdifft2 = work.kdifft2;
% old_diff.kdiffq1 = work.kdiffq1;
% old_diff.kdiffq2 = work.kdiffq2;
% old_diff.kdiffu1 = work.kdiffu1;
% old_diff.kdiffu2 = work.kdiffu2;
% old_diff.kdiffw1 = work.kdiffw1;
% old_diff.kdiffw2 = work.kdiffw2;
% old_diff.kdifftke1 = work.kdifftke1;
% old_diff.kdifftke2 = work.kdifftke2;


% Right hand sides of all equations
bdt = time.beta*time.dt;
if switches.c
    rhs1m   = state_old.fluid(2).m          + state_old.fluid(1).m          + bdt*tend.fluid(1).m.tot;
    rhs1eta = state_old.fluid(2).eta.*m2bar + state_old.fluid(1).eta.*m1bar + bdt*tend.fluid(1).meta.tot;
    rhs1q   = state_old.fluid(2).q  .*m2bar + state_old.fluid(1).q  .*m1bar + bdt*tend.fluid(1).mq.tot;
    rhs1w   = state_old.fluid(2).w  .*m2bar + state_old.fluid(1).w  .*m1bar + bdt*tend.fluid(1).mw.tot;
    rhs1u   = state_old.fluid(2).u  .*m2    + state_old.fluid(1).u  .*m1    + bdt*tend.fluid(1).mu.tot;
    rhs1v   = state_old.fluid(2).v  .*m2    + state_old.fluid(1).v  .*m1    + bdt*tend.fluid(1).mv.tot;
    rhs1tke = state_old.fluid(2).tke.*m2    + state_old.fluid(1).tke.*m1    + bdt*tend.fluid(1).mtke.tot;
    rhs2m   = bdt*tend.fluid(2).m.tot;
    rhs2eta = bdt*tend.fluid(2).meta.tot;
    rhs2q   = bdt*tend.fluid(2).mq.tot;
    rhs2w   = bdt*tend.fluid(2).mw.tot;
    rhs2u   = bdt*tend.fluid(2).mu.tot;
    rhs2v   = bdt*tend.fluid(2).mv.tot;
    rhs2tke = bdt*tend.fluid(2).mtke.tot;
else
    rhs1m   = state_old.fluid(1).m          + bdt*tend.fluid(1).m.tot;
    rhs1eta = state_old.fluid(1).eta.*m1bar + bdt*tend.fluid(1).meta.tot;
    rhs1q   = state_old.fluid(1).q  .*m1bar + bdt*tend.fluid(1).mq.tot;
    rhs1w   = state_old.fluid(1).w  .*m1bar + bdt*tend.fluid(1).mw.tot;
    rhs1u   = state_old.fluid(1).u  .*m1    + bdt*tend.fluid(1).mu.tot;
    rhs1v   = state_old.fluid(1).v  .*m1    + bdt*tend.fluid(1).mv.tot;
    rhs1tke = state_old.fluid(1).tke.*m1    + bdt*tend.fluid(1).mtke.tot;
    rhs2m   = state_old.fluid(2).m          + bdt*tend.fluid(2).m.tot;
    rhs2eta = state_old.fluid(2).eta.*m2bar + bdt*tend.fluid(2).meta.tot;
    rhs2q   = state_old.fluid(2).q  .*m2bar + bdt*tend.fluid(2).mq.tot;
    rhs2w   = state_old.fluid(2).w  .*m2bar + bdt*tend.fluid(2).mw.tot;
    rhs2u   = state_old.fluid(2).u  .*m2    + bdt*tend.fluid(2).mu.tot;
    rhs2v   = state_old.fluid(2).v  .*m2    + bdt*tend.fluid(2).mv.tot;
    rhs2tke = state_old.fluid(2).tke.*m2    + bdt*tend.fluid(2).mtke.tot;
end

% Accumulate forcing and budget terms
accdt = bdt;
accumulate

% ------

if istep == 10000000
    
    % Do a convergence test
    % First brute force a converged solution
    qn_iter_max = 20;
    conv_diag = 0;
    quasinewton
    
    % Switch on approximation to be tested
    %switches.c = 1
    % and diagnose residual fixer to mimic converged solution
    find_res_fix
    
    % Save converged solution
    state_tru = state_new;
    tend_tru  = tend;
    eos_tru = eos;
    work_tru = work;
    
    % Reset and perturb first guess
    perturb_state
    disp('First guess has been reset')
    pause
    
    % Now go again with convergence diagnostics
    qn_iter_max = 10;
    conv_diag = 1;
    quasinewton
    pause
    
else
    
    % Standard time step
    % Number of quasi-Newton iterations
    qn_iter_max = 4;

    % Carry out quasi-Newton iterations
    conv_diag = 0;
    quasinewton

end

% Fix descent in updrafts by homogenizing fluids 1 and 2
[inc_w1,inc_w2,inc_eta1,inc_eta2,inc_q1,inc_q2] = fix_negative_w(grid,state_new);
% disp('*** no w fixer ***')
state_new.fluid(1).w   = state_new.fluid(1).w   + inc_w1;
state_new.fluid(2).w   = state_new.fluid(2).w   + inc_w2;
state_new.fluid(1).eta = state_new.fluid(1).eta + inc_eta1;
state_new.fluid(2).eta = state_new.fluid(2).eta + inc_eta2;
state_new.fluid(1).q   = state_new.fluid(1).q   + inc_q1;
state_new.fluid(2).q   = state_new.fluid(2).q   + inc_q2;

% and save increments for budgets
budgets.w2.wfix   = inc_w2/time.dt;
budgets.eta1.wfix = inc_eta1./time.dt;
budgets.eta2.wfix = inc_eta2./time.dt;
budgets.q1.wfix   = inc_q1./time.dt;
budgets.q2.wfix   = inc_q2./time.dt;


% Update old time level values ready for next step
state_old = state_new;
