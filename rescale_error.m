% Rescale the errors within the quasi-Newton solver

% First compute the errors
state_err.p            = state_new.p            - state_tru.p;
state_err.fluid(1).m   = state_new.fluid(1).m   - state_tru.fluid(1).m;
state_err.fluid(1).w   = state_new.fluid(1).w   - state_tru.fluid(1).w;
state_err.fluid(1).u   = state_new.fluid(1).u   - state_tru.fluid(1).u;
state_err.fluid(1).v   = state_new.fluid(1).v   - state_tru.fluid(1).v;
state_err.fluid(1).tke = state_new.fluid(1).tke - state_tru.fluid(1).tke;
state_err.fluid(1).eta = state_new.fluid(1).eta - state_tru.fluid(1).eta;
state_err.fluid(1).q   = state_new.fluid(1).q   - state_tru.fluid(1).q;
state_err.fluid(1).T   = state_new.fluid(1).T   - state_tru.fluid(1).T;
state_err.fluid(1).Tw  = state_new.fluid(1).Tw  - state_tru.fluid(1).Tw;
state_err.fluid(2).m   = state_new.fluid(2).m   - state_tru.fluid(2).m;
state_err.fluid(2).w   = state_new.fluid(2).w   - state_tru.fluid(2).w;
state_err.fluid(2).u   = state_new.fluid(2).u   - state_tru.fluid(2).u;
state_err.fluid(2).v   = state_new.fluid(2).v   - state_tru.fluid(2).v;
state_err.fluid(2).tke = state_new.fluid(2).tke - state_tru.fluid(2).tke;
state_err.fluid(2).eta = state_new.fluid(2).eta - state_tru.fluid(2).eta;
state_err.fluid(2).q   = state_new.fluid(2).q   - state_tru.fluid(2).q;
state_err.fluid(2).T   = state_new.fluid(2).T   - state_tru.fluid(2).T;
state_err.fluid(2).Tw  = state_new.fluid(2).Tw  - state_tru.fluid(2).Tw;

% disp('Errors before rescaling')
% disp(['err_p     ' num2str(state_err.p(krange))])
% disp(['err_w1    ' num2str(state_err.fluid(1).w(krange))])
% disp(['err_w2    ' num2str(state_err.fluid(2).w(krange))])
% disp(['err_m1    ' num2str(state_err.fluid(1).m(krange))])
% disp(['err_m2    ' num2str(state_err.fluid(2).m(krange))])
% disp(['err_eta1  ' num2str(state_err.fluid(1).eta(krange))])
% disp(['err_eta2  ' num2str(state_err.fluid(2).eta(krange))])
% disp(['err_q1    ' num2str(state_err.fluid(1).q(krange))])
% disp(['err_q2    ' num2str(state_err.fluid(2).q(krange))])
% disp(['err_T1    ' num2str(state_err.fluid(1).T(krange))])
% disp(['err_T2    ' num2str(state_err.fluid(2).T(krange))])
% disp(['err_Tw1   ' num2str(state_err.fluid(1).Tw(krange))])
% disp(['err_Tw2   ' num2str(state_err.fluid(2).Tw(krange))])


% Determine a measure of error amplitude
err_amp = max(abs(state_err.p)) + 1e-12;
disp(['original p error ' num2str(state_err.p(krange))])

% Rescaling factor
%f_rescale = 0.1/err_amp
f_rescale = 1.0
%f_rescale0 = f_rescale;
f_rescale0 = 0.0;
if qn_iter == 12
    f_rescale1 = 0.0;
else
    f_rescale1 = f_rescale;
end


% Rescale errors
state_err.p            = f_rescale *state_err.p;
state_err.fluid(1).m   = f_rescale *state_err.fluid(1).m;
state_err.fluid(1).w   = f_rescale *state_err.fluid(1).w;
state_err.fluid(1).u   = f_rescale *state_err.fluid(1).u;
state_err.fluid(1).v   = f_rescale *state_err.fluid(1).v;
state_err.fluid(1).tke = f_rescale *state_err.fluid(1).tke;
state_err.fluid(1).eta = f_rescale *state_err.fluid(1).eta;
state_err.fluid(1).q   = f_rescale *state_err.fluid(1).q;
state_err.fluid(1).T   = f_rescale *state_err.fluid(1).T;
state_err.fluid(1).Tw  = f_rescale *state_err.fluid(1).Tw;
state_err.fluid(2).m   = f_rescale *state_err.fluid(2).m;
state_err.fluid(2).w   = f_rescale *state_err.fluid(2).w;
state_err.fluid(2).u   = f_rescale *state_err.fluid(2).u;
state_err.fluid(2).v   = f_rescale *state_err.fluid(2).v;
state_err.fluid(2).tke = f_rescale *state_err.fluid(2).tke;
state_err.fluid(2).eta = f_rescale *state_err.fluid(2).eta;
state_err.fluid(2).q   = f_rescale *state_err.fluid(2).q;
state_err.fluid(2).T   = f_rescale *state_err.fluid(2).T;
state_err.fluid(2).Tw  = f_rescale *state_err.fluid(2).Tw;

disp('Errors after rescaling')
disp(['err_p     ' num2str(state_err.p(krange))])
disp(['err_w1    ' num2str(state_err.fluid(1).w(krange))])
disp(['err_w2    ' num2str(state_err.fluid(2).w(krange))])
disp(['err_u1    ' num2str(state_err.fluid(1).u(krange))])
disp(['err_u2    ' num2str(state_err.fluid(2).u(krange))])
disp(['err_v1    ' num2str(state_err.fluid(1).v(krange))])
disp(['err_v2    ' num2str(state_err.fluid(2).v(krange))])
disp(['err_m1    ' num2str(state_err.fluid(1).m(krange))])
disp(['err_m2    ' num2str(state_err.fluid(2).m(krange))])
disp(['err_eta1  ' num2str(state_err.fluid(1).eta(krange))])
disp(['err_eta2  ' num2str(state_err.fluid(2).eta(krange))])
disp(['err_q1    ' num2str(state_err.fluid(1).q(krange))])
disp(['err_q2    ' num2str(state_err.fluid(2).q(krange))])
disp(['err_tke1  ' num2str(state_err.fluid(1).tke(krange))])
disp(['err_tke2  ' num2str(state_err.fluid(2).tke(krange))])
disp(['err_T1    ' num2str(state_err.fluid(1).T(krange))])
disp(['err_T2    ' num2str(state_err.fluid(2).T(krange))])
disp(['err_Tw1   ' num2str(state_err.fluid(1).Tw(krange))])
disp(['err_Tw2   ' num2str(state_err.fluid(2).Tw(krange))])

% And update state
state_new.p            = state_tru.p            + state_err.p;
state_new.fluid(1).m   = state_tru.fluid(1).m   + state_err.fluid(1).m;
state_new.fluid(1).w   = state_tru.fluid(1).w   + state_err.fluid(1).w;
state_new.fluid(1).u   = state_tru.fluid(1).u   + state_err.fluid(1).u;
state_new.fluid(1).v   = state_tru.fluid(1).v   + state_err.fluid(1).v;
state_new.fluid(1).tke = state_tru.fluid(1).tke + state_err.fluid(1).tke;
state_new.fluid(1).eta = state_tru.fluid(1).eta + state_err.fluid(1).eta;
state_new.fluid(1).q   = state_tru.fluid(1).q   + state_err.fluid(1).q;
state_new.fluid(1).T   = state_tru.fluid(1).T   + state_err.fluid(1).T;
state_new.fluid(1).Tw  = state_tru.fluid(1).Tw  + state_err.fluid(1).Tw;
state_new.fluid(2).m   = state_tru.fluid(2).m   + state_err.fluid(2).m;
state_new.fluid(2).w   = state_tru.fluid(2).w   + state_err.fluid(2).w;
state_new.fluid(2).u   = state_tru.fluid(2).u   + state_err.fluid(2).u;
state_new.fluid(2).v   = state_tru.fluid(2).v   + state_err.fluid(2).v;
state_new.fluid(2).tke = state_tru.fluid(2).tke + state_err.fluid(2).tke;
state_new.fluid(2).eta = state_tru.fluid(2).eta + state_err.fluid(2).eta;
state_new.fluid(2).q   = state_tru.fluid(2).q   + state_err.fluid(2).q;
state_new.fluid(2).T   = state_tru.fluid(2).T   + state_err.fluid(2).T;
state_new.fluid(2).Tw  = state_tru.fluid(2).Tw  + state_err.fluid(2).Tw;


figure(20)
subplot(2,3,1)
plot(state_err.p(prange),grid.zp(prange))
title(['p err  f =',num2str(f_rescale)])
set(gca,'fontsize',fs)
subplot(2,3,2)
plot(state_err.fluid(1).w(prange),grid.zw(prange),'b',state_err.fluid(2).w(prange),grid.zw(prange),'r')
title('w err')
set(gca,'fontsize',fs)
subplot(2,3,3)
plot(state_err.fluid(1).m(prange),grid.zp(prange),'b',state_err.fluid(2).m(prange),grid.zp(prange),'r')
title('m err')
set(gca,'fontsize',fs)
subplot(2,3,4)
plot(state_err.fluid(1).eta(prange),grid.zw(prange),'b',state_err.fluid(2).eta(prange),grid.zw(prange),'r')
title('eta err')
set(gca,'fontsize',fs)
subplot(2,3,5)
plot(state_err.fluid(1).q(prange),grid.zw(prange),'b',state_err.fluid(2).q(prange),grid.zw(prange),'r')
title('q err')
set(gca,'fontsize',fs)
subplot(2,3,6)
plot(state_err.fluid(1).T(prange),grid.zp(prange),'b',state_err.fluid(2).T(prange),grid.zp(prange),'r',...
state_err.fluid(1).Tw(prange),grid.zw(prange),'b--',state_err.fluid(2).Tw(prange),grid.zw(prange),'r--')
title('T err')
set(gca,'fontsize',fs)
% pause %(0.8)
figure(1)


