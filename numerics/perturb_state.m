% Perturb first guess for new time level state in order to
% perform a convergence test

% ------

% Either perturb a specific field...

kpert = 2;
%kpert = 1:nz;
%disp(['p pertn  '  num2str(0.001*state_new.p(kpert))])
%state_new.p(kpert) = state_new.p(kpert)*1.001;
%disp(['p pertn  '  num2str(0.01)])
%state_new.p(kpert) = state_new.p(kpert) + 0.01;
%disp(['m1 pertn  '  num2str(0.001)])
%state_new.fluid(1).m(kpert) = state_new.fluid(1).m(kpert) + 0.001;
%disp(['m1 pertn  '  num2str(0.01*state_new.fluid(1).m(kpert))])
%state_new.fluid(1).m(kpert) = state_new.fluid(1).m(kpert)*1.01;
%disp(['m2 pertn  '  num2str(0.01*state_new.fluid(2).m(kpert))])
%disp(['m2 pertn  '  num2str(0.001)])
%state_new.fluid(2).m(kpert) = state_new.fluid(2).m(kpert) - 0.001;
%disp(['m2 pertn  '  num2str( 0.01)])
%state_new.fluid(2).m(kpert) = state_new.fluid(2).m(kpert) + 0.01;
%disp(['m1 m2 pertn  '  num2str(0.001)])
%state_new.fluid(1).m(kpert) = state_new.fluid(1).m(kpert) - 0.001;
%state_new.fluid(2).m(kpert) = state_new.fluid(2).m(kpert) - 0.01;
%state_new.fluid(1).m        = state_new.fluid(1).m        - 0.001;
%state_new.fluid(2).m        = state_new.fluid(2).m        - 0.01;
%disp(['w1 pertn  '  num2str(0.001)])
%state_new.fluid(1).w(kpert) = state_new.fluid(1).w(kpert) - 0.001;
%disp(['w2 pertn  '  num2str(0.001)])
%state_new.fluid(2).w(kpert) = state_new.fluid(2).w(kpert) - 0.001;
%disp(['eta1 pertn  '  num2str(0.01*state_new.fluid(1).eta(kpert))])
%state_new.fluid(1).eta(kpert) = state_new.fluid(1).eta(kpert)*1.01;
%disp(['eta2 pertn  '  num2str(0.01*state_new.fluid(2).eta(kpert))])
%state_new.fluid(2).eta(kpert) = state_new.fluid(2).eta(kpert)*0.99;
%disp(['q1 pertn  '  num2str(0.0001)])
%state_new.fluid(1).q(kpert) = state_new.fluid(1).q(kpert) + 0.0001;
%disp(['q2 pertn  '  num2str(0.0001)])
%state_new.fluid(2).q(kpert) = state_new.fluid(2).q(kpert) - 0.0001;
%disp(['T1 pertn  '  num2str(0.0005*state_new.fluid(1).T(kpert))])
%state_new.fluid(1).T(kpert) = state_new.fluid(1).T(kpert)*1.0005;
%disp(['T2 pertn  '  num2str(0.002*state_new.fluid(2).T(kpert))])
%state_new.fluid(2).T(kpert) = state_new.fluid(2).T(kpert)*1.002;
%disp(['Tw1 pertn  '  num2str(0.0005*state_new.fluid(1).Tw(kpert))])
%state_new.fluid(1).Tw(kpert) = state_new.fluid(1).Tw(kpert)*1.0005;
%disp(['Tw2 pertn  '  num2str(0.002*state_new.fluid(2).Tw(kpert))])
%state_new.fluid(2).Tw(kpert) = state_new.fluid(2).Tw(kpert)*1.002;
%disp(['u1 pertn  '  num2str(0.1)])
%state_new.fluid(1).u(kpert) = state_new.fluid(1).u(kpert) + 0.1;
%disp(['u2 pertn  '  num2str(0.1)])
%state_new.fluid(2).u(kpert) = state_new.fluid(2).u(kpert) + 0.1;

% ------

% Or reset to previous time level

disp('Reset to old state')
state_new = state_old;

% ------

% Or leave converged solution - do nothing!

% ------
