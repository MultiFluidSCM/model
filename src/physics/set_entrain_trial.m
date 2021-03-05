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
% 4     Instability source + combination of mixing and sorting detrainment
%                            sorting determines amount and properties
%                          + dynamic entrainment/detrainment based on dw/dz
ischeme = 4;

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
tke1 = state.fluid(1).tke;
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
dw1dz = (w1(2:nzp) - w1(1:nz))./grid.dzp;
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


% Instability source entrainment where stratification is unstable
set_entrain_trial_instab

% Turbulent mixing entrainment & detrainment
set_entrain_trial_mix

% APDF-based formulation for detrainment, sorting pdf
set_entrain_trial_sort

% Dynamical entrainment/detrainment based on vertical velocity convergence
set_entrain_trial_dwdz

% Remove new transfer contributions for old schemes
if ischeme == 0
    relabel.M12_instab = 0*relabel.M12_instab;
    relabel.M12_sort = 0*relabel.M12_sort;
elseif ischeme == 1 | ischeme == 3
    relabel.M12_instab = 0*relabel.M12_instab;
end

% Combined entrainment (M21) and detrainment (M12) rates
relabel.M12 = relabel.M12_instab + relabel.M12_mix + relabel.M12_sort + relabel.M12_dwdz;
relabel.M21 = relabel.M21_instab + relabel.M21_mix + relabel.M21_sort + relabel.M21_dwdz;

% Get the properties of the entrained and detrained fluids
set_entrain_trial_transfer_properties


if ischeme == 0 | ischeme == 1
    % Safe option for detrained values
    relabel.etahat12 = relabel.etahat12_mix;
    relabel.qhat12   = relabel.qhat12_mix;
    relabel.what12   = relabel.what12_mix;
    % Factor needed for linearized variance equation
    relabel.f_sort_chi_hat = 0*chi_hat;
elseif ischeme == 3 | ischeme == 4
    denominator = relabel.M12 + 1e-8*(relabel.M12 == 0);
    f_sort = weight_to_w(grid,relabel.M12_sort./denominator);
    
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
relabel.dM21dm1   = dM21dm1_mix + dM21dm1_instab + dM21dm1_dwdz;% + dM21dm1_sort; % Is sort contribution already included in Newton solver separately?
relabel.dM21dm2   = dM21dm2_mix + dM21dm2_instab + dM21dm2_dwdz;% + dM21dm2_sort;
relabel.dM21dw1   = zeros(1,nz);
relabel.dM21dw2   = zeros(1,nz);
relabel.dM21deta1 = zeros(1,nz);
relabel.dM21deta2 = zeros(1,nz);
relabel.dM21dq1   = zeros(1,nz);
relabel.dM21dq2   = zeros(1,nz);

relabel.dM12dm1   = dM12dm1_mix + dM12dm1_instab + dM12dm1_dwdz;% + dM12dm1_sort;
relabel.dM12dm2   = dM12dm2_mix + dM12dm2_instab + dM12dm2_dwdz;% + dM12dm2_sort;
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

% Save background profile 
relabel.ideal = ones(1,nz)*sigma20;

end