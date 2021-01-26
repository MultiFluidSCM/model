function relabel = set_entrain_trial(grid, state, buoy, eos, scales, kdiffw2, constants, dt)

% Set mass entrainment and detrainment rates
% M_ij = ( m_i m_j / rho_i ) * rate_ij
%
% the values of all entrained fields what, etahat, qhat, uhat, vhat
%
% and their derivatives

% For reproducibility, define certain `standard' schemes
% 0     Instability source + relaxation to desired profile
% 1     Instability source + combination of mixing and sorting detrainment
%                            sorting determines amount but not properties
% 3     Instability source + combination of mixing and sorting detrainment
%                            sorting determines amount and properties
ischeme = 3;

% First unpack the fields needed
nz = grid.nz;
nzp = nz + 1;
zp = grid.zp;
abovep = grid.abovep;
belowp = grid.belowp;
m1   = state.fluid(1).m;
m2   = state.fluid(2).m;
eta1 = state.fluid(1).eta;
eta2 = state.fluid(2).eta;
w1   = state.fluid(1).w;
w2   = state.fluid(2).w;
q1   = state.fluid(1).q;
q2   = state.fluid(2).q;
u1   = state.fluid(1).u;
u2   = state.fluid(2).u;
v1   = state.fluid(1).v;
v2   = state.fluid(2).v;
tke2 = state.fluid(2).tke;
bentraint = constants.param.bentraint;
bentrainq = constants.param.bentrainq;
bentrainw = constants.param.bentrainw;
bentrainu = constants.param.bentrainu;
bdetraint = constants.param.bdetraint;
bdetrainq = constants.param.bdetrainq;
bdetrainw = constants.param.bdetrainw;
bdetrainu = constants.param.bdetrainu;
sort = constants.param.sort;
dwdz = constants.param.dwdz;
mix = constants.param.mix;
instab = constants.param.instab;

sigma1 = eos.sigma1;
sigma2 = eos.sigma2;

% ------

% Preliminary calculations

% Background sigma2 when nothing is happening
sigma20 = constants.param.sigma00;
sigma10 = 1 - sigma20;

% To bound relaxation rate
rdt = 1/dt;

% Buoyancy on p levels
buoybar = (abovep.*buoy(2:nzp) + belowp.*buoy(1:nz));

% Mean density
rho = m1 + m2;

% Mass fractions
mf1 = m1./rho;
mf2 = 1 - mf1;

% Inverse sqrt of 2
rroot2 = sqrt(0.5);

% sqrt of -1 * buoyancy frequency squared (where it's unstable)
% or buoyancy frequency (where its stable)
n1neg = min(sqrt(max(-eos.nsq1,0)),rdt);
n1pos = min(sqrt(max( eos.nsq1,0)),rdt);
n2neg = min(sqrt(max(-eos.nsq2,0)),rdt);
n2pos = min(sqrt(max( eos.nsq2,0)),rdt);

% Vertically averaged updraft w
w2bar = max(0,grid.abovep.*w2(2:nzp) + grid.belowp.*w2(1:nz));
% and vertical derivative
dw2dz = (w2(2:nzp) - w2(1:nz))./grid.dzp;

% Subfilter standard deviation of w, assuming tke has equal contributions
% from u, v and w
wstd = sqrt(tke2*2/3);
wstdw(2:nz) = grid.abovew(2:nz).*wstd(2:nz) + grid.beloww(2:nz).*wstd(1:nz-1);
wstdw(1) = wstdw(2);
wstdw(nzp) = wstdw(nz);
% And subfilter standard deviations of eta, q, and buoyancy
etastd = sqrt(state.fluid(2).vareta);
qstd   = sqrt(state.fluid(2).varq);
bstd = constants.phys.gravity*sqrt(eos.Varrhow2)./eos.rhow2;


% shear between fluids
du = u2 - u1;
dv = v2 - v1;
dww = w2 - w1;
dw = grid.abovep.*dww(2:nzp) + grid.belowp.*dww(1:nz);
absdu = sqrt(du.*du + dv.*dv + dw.*dw);

% ------


% Instability source entrainment where stratification is unstable
% *** Put tuneable coefficients in constants.params ***
relabel.M21_instab = 0.2*m1.*n1neg;
dM21dm1_instab = 0.2*n1neg;


% Relaxation to homogeneous state with background sigma
% Ratio of sigma2 to reference value
rr = sigma2/sigma20;
% Merging factor to go smoothly between fast relaxation rate r1 for rr < 1
% to turbulent rate r2 for rr > 2
rr = min(2, max(1,rr) );
qq = (rr - 2).*(rr - 2).*(1 - (rr - 1).*(rr - 1));
r1 = rdt;
if ischeme == 0
    r2 = n1pos;
elseif ischeme == 1
    r2 = min(0.5*sqrt(tke2)./scales.L_plume,rdt);
elseif ischeme == 3
    r2 = min(0.25*sqrt(tke2)./scales.L_plume,rdt);
else
    disp('unknown scheme in set_entrain_trial')
    pause
end
rate_mix = r2 + (r1 - r2).*qq;
relabel.rate_mix = rate_mix; % Just for diagnostics
% Reference profiles of  m1 and m2
m2ref = sigma20*rho;
m1ref = rho - m2ref;
% Factors needed to compute mixing
m1fac = 2*mf1 - sigma10;
m2fac = 2*mf2 - sigma20;
% There are three cases:
% m1fac < 0, 0 < m1fac < 1, 1 < m1fac
case1 = m1fac < 0;
case3 = 1 < m1fac;
case2 = 1 - (case1 + case3);
% Bound to prevent negative relabelling rates
m1fac = min(1, max(0,m1fac) );
m2fac = min(1, max(0,m2fac) );
% Relabelling rates
relabel.M21_mix = rate_mix.*m2.*m1fac;
relabel.M12_mix = rate_mix.*m1.*m2fac;
% Derivatives *** Need to think carefully about these ***
dM21dm1_mix = case1.*zeros(1,nz) + case2.*2.*rate_mix.*mf2.*mf2            + case3.*zeros(1,nz);
dM21dm2_mix = case1.*zeros(1,nz) + case2.*rate_mix.*(2*mf1.*mf1 - sigma10) + case3.*rate_mix;
dM12dm1_mix = case1.*rate_mix    + case2.*rate_mix.*(2*mf2.*mf2 - sigma20) + case3.*zeros(1,nz);
dM12dm2_mix = case1.*zeros(1,nz) + case2.*2.*rate_mix.*mf1.*mf1            + case3.*zeros(1,nz);




% APDF-based formulation for detrainment, sorting pdf
relabel.w2bar = w2bar;
relabel.wstd = wstd;
% 1/sqrt(2) * normalized variable at cut of pdf on p-levels
chi_cut = - rroot2*w2bar./wstd;
% and w-levels (bounded to prevent divide by zero below)
% chi_cut_w = max(-rroot2*w2./wstdw,-4);
chi_cut_w = -rroot2*w2./wstdw;
% Fraction of pdf cut off i.e. with w < 0
relabel.frac = 0.5*(1 + erf(chi_cut));
% disp(['dw2dz = ',num2str(dw2dz(52))])
% if abs(dw2dz(52) + 0.0023576) < 1e-6
%     save('dw2dz.mat','dw2dz')
%     disp('dw2dz saved')
% else
%     load('dw2dz.mat')
%     disp('dw2dz loaded')
% end
% Rate at which it is cut off
% rate_sort = 5*(wstd./scales.L_plume);
if ischeme == 1
    rate_sort = min(10*max(0,-dw2dz),rdt);
elseif ischeme == 3
    %rate_sort = min(20*max(0,-dw2dz),rdt);
    disp('*** Experimental rate_sort, half rate ***')
    % rate_sort = sign(n1pos).*min(20*rate_mix,rdt);  % Strictly this is scheme 7 or 8
    % rate_sort = (buoybar < 0).*min(20*rate_mix,rdt);  % 
    % ss3 = 0.5*(1 - tanh( 2*buoy./(bstd + eps) + 0.4));  % Smooth switch 10b
    % ss3 = 0.5*(1 - tanh( 2*buoy./(bstd + eps) + 0.8));  % Smooth switch 10d
    % ss3 = 0.5*(1 - tanh( 5*buoy./(bstd + eps) + 2));  % Smooth switch 10c
    bnorm = weight_to_w(grid,tke2./scales.L_turb2);
    ss3 = 0.5*(1 - tanh( 3*buoy./bnorm + 1.2));  % Smooth switch 10e
    relabel.ss3bar = abovep.*ss3(2:nzp) + belowp.*ss3(1:nz);
    rate_sort = relabel.ss3bar.*min(10*rate_mix,rdt);  % 
    % rate_buoy = max(-buoybar,0)./wstd;
    % rate_sort = min(20*rate_buoy,rdt);
else
    rate_sort = 0;
end
% Resulting detrainment rate
relabel.M12_sort = rate_sort.*m2.*relabel.frac;
% disp(['wstd         = ',num2str(wstd(69))])
% disp(['buoybar      = ',num2str(buoybar(69))])
% disp(['rate_sort    = ',num2str(rate_sort(69))])
% disp(['rate_mix     = ',num2str(rate_mix(69))])
% disp(['m2           = ',num2str(m2(69))])
% disp(['relabel.frac = ',num2str(relabel.frac(69))])


% w, eta and q of detrained fluid
% chi_hat = - sqrt(2/pi)*exp(-chi_cut_w.^2)./(1 + erf(chi_cut_w));
for k = 1:nzp
    [chi_hat(k),deriv(k)] = compute_chi_hat(chi_cut_w(k));
end
% if abs(chi_hat(52) + 0.896576265536613) < 1e-8
%     tempeta2 = eta2;
%     save('chihat.mat','chi_hat','tempeta2','etastd')
%     disp('chi_hat saved')
% else
%     load('chihat.mat')
%     disp('chi_hat loaded')
% end
relabel.what12_sort   = w2   + wstdw .*chi_hat;
relabel.etahat12_sort = eta2 + etastd.*chi_hat;
%relabel.etahat12_sort = tempeta2 + etastd.*chi_hat;    % ****
relabel.qhat12_sort   = q2   + qstd  .*chi_hat;
% Derivative of what12_sort wrt wstdw
relabel.what_deriv = deriv;



% Entrained and detrained values of eta
[relabel.etahat12_mix,relabel.detahat12deta1,relabel.detahat12deta2,...
 relabel.etahat21    ,relabel.detahat21deta1,relabel.detahat21deta2] ...
    = findqhat(eta1,eta2,bentraint,bdetraint);

% Entrained and detrained values of q
[relabel.qhat12_mix,relabel.dqhat12dq1,relabel.dqhat12dq2,...
 relabel.qhat21    ,relabel.dqhat21dq1,relabel.dqhat21dq2] ...
    = findqhat(q1,q2,bentrainq,bdetrainq);

% Entrained and detrained values of w
[relabel.what12_mix,relabel.dwhat12dw1,relabel.dwhat12dw2,...
 relabel.what21_mix,relabel.dwhat21dw1,relabel.dwhat21dw2] ...
    = findqhat(w1,w2,mix.bentrainw,mix.bdetrainw);
relabel.what21_mix = max(relabel.what21_mix,0);

% Entrained and detrained values of u and v
[relabel.uhat12,relabel.duhat12du1,relabel.duhat12du2,...
 relabel.uhat21,relabel.duhat21du1,relabel.duhat21du2] ...
    = findqhat(u1,u2,bentrainu,bdetrainu);
[relabel.vhat12,relabel.dvhat12dv1,relabel.dvhat12dv2,...
 relabel.vhat21,relabel.dvhat21dv1,relabel.dvhat21dv2] ...
    = findqhat(v1,v2,bentrainu,bdetrainu);
    

% Entrained and detrained values of w
[relabel.what12_instab,relabel.dwhat12dw1,relabel.dwhat12dw2,...
 relabel.what21_instab,relabel.dwhat21dw1,relabel.dwhat21dw2] ...
    = findqhat(w1,w2,instab.bentrainw,instab.bdetrainw);
relabel.what21_instab = max(relabel.what21_instab,0);
[relabel.what12_dwdz,relabel.dwhat12dw1,relabel.dwhat12dw2,...
 relabel.what21_dwdz,relabel.dwhat21dw1,relabel.dwhat21dw2] ...
    = findqhat(w1,w2,dwdz.bentrainw,dwdz.bdetrainw);
relabel.what21_dwdz = max(relabel.what21_dwdz,0);

% Combined entrainment rate
relabel.M21 = relabel.M21_instab + relabel.M21_mix;

denominator = relabel.M21 + 1e-8*(relabel.M21 == 0);
f_instab = weight_to_w(grid,relabel.M21_instab./denominator);
f_mix  = 1 - f_instab;
relabel.what12_blend   = f_instab.*relabel.what12_instab ...
                       + f_mix .*relabel.what12_mix;
relabel.what12 = relabel.what12_blend;


% For testing
if ischeme == 0
    relabel.M12 = relabel.M12_mix;
elseif ischeme == 1 | ischeme == 3
    relabel.M12 = relabel.M12_sort + relabel.M12_mix;
else
    disp('unknown scheme in set_entrain_trial')
    pause
end


% Experimental option for detrained values
% To catch divide by zero
denominator = relabel.M12 + 1e-8*(relabel.M12 == 0);
f_sort = weight_to_w(grid,relabel.M12_sort./denominator);
f_mix  = 1 - f_sort;
relabel.f_sort = f_sort;
relabel.etahat12_blend = f_sort.*relabel.etahat12_sort ...
                       + f_mix .*relabel.etahat12_mix;
relabel.qhat12_blend   = f_sort.*relabel.qhat12_sort ...
                       + f_mix .*relabel.qhat12_mix;
relabel.what12_blend   = f_sort.*relabel.what12_sort ...
                       + f_mix .*relabel.what12_mix;

                   
if ischeme == 0 | ischeme == 1
    % Safe option for detrained values
    relabel.etahat12 = relabel.etahat12_mix;
    relabel.qhat12   = relabel.qhat12_mix;
    relabel.what12   = relabel.what12_mix;
    % Factor needed for linearized variance equation
    relabel.f_sort_chi_hat = 0*chi_hat;
elseif ischeme == 3
    % Apply experimental option
    relabel.etahat12 = relabel.etahat12_blend;
    relabel.qhat12   = relabel.qhat12_blend;
    relabel.what12   = relabel.what12_blend;
    % Factor needed for linearized variance equation
    relabel.f_sort_chi_hat = f_sort.*chi_hat;
else
    disp('unknown scheme in set_entrain_trial')
    pause
end

% Make sure what's are zero at top and bottom boundaries
relabel.what12(1  )   = 0;
relabel.what12(nzp)   = 0;
relabel.what21(1  )   = 0;
relabel.what21(nzp)   = 0;
    
% Derivatives
relabel.dM21dm1   = dM21dm1_mix + dM21dm1_instab;
relabel.dM21dm2   = dM21dm2_mix;
relabel.dM21dw1   = zeros(1,nz);
relabel.dM21dw2   = zeros(1,nz);
relabel.dM21deta1 = zeros(1,nz);
relabel.dM21deta2 = zeros(1,nz);
relabel.dM21dq1   = zeros(1,nz);
relabel.dM21dq2   = zeros(1,nz);

relabel.dM12dm1   = dM12dm1_mix;
relabel.dM12dm2   = dM12dm2_mix;
relabel.dM12dw1   = zeros(1,nz);
relabel.dM12dw2   = zeros(1,nz);
relabel.dM12deta1 = zeros(1,nz);
relabel.dM12deta2 = zeros(1,nz);
relabel.dM12dq1   = zeros(1,nz);
relabel.dM12dq2   = zeros(1,nz);
    


% Interpolate entrainment and detrainment to w levels
% using `reversed' weighting for conservation
relabel.M12bar = weight_to_w(grid,relabel.M12);
relabel.M21bar = weight_to_w(grid,relabel.M21);

% Save bckground profile 
relabel.ideal = ones(1,nz)*sigma20;

% relabel.M12 = M12;
% relabel.M21 = M21;
% relabel.M12bar = M12bar;
% relabel.M21bar = M21bar;
% relabel.dM21dm1 = dM21dm1;
% relabel.dM21dm2 = dM21dm2;
% relabel.dM21dw1 = dM21dw1;
% relabel.dM21dw2 = dM21dw2;
% relabel.dM21deta1 = dM21deta1;
% relabel.dM21deta2 = dM21deta2;
% relabel.dM21dq1 = dM21dq1;
% relabel.dM21dq2 = dM21dq2;
% relabel.dM12dm1 = dM12dm1;
% relabel.dM12dm2 = dM12dm2;
% relabel.dM12dw1 = dM12dw1;
% relabel.dM12dw2 = dM12dw2;
% relabel.dM12deta1 = dM12deta1;
% relabel.dM12deta2 = dM12deta2;
% relabel.dM12dq1 = dM12dq1;
% relabel.dM12dq2 = dM12dq2;




end