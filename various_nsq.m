% Various ways of computing buoyancy frequency squared


% Vertical pressure gradient
dpdz(2:nz)   = (state.p(2:nz) - state.p(1:nz-1))./grid.dzw(2:nz);
dpdz(1) = dpdz(2);
dpdz(nzp) = dpdz(nz);
dpdzbar = grid.abovep.*dpdz(2:nzp) + grid.belowp.*dpdz(1:nz);

% Vertical eta and q gradients
deta1dz = (state.fluid(1).eta(2:nzp) - state.fluid(1).eta(1:nz))./grid.dzp;
dq1dz = (state.fluid(1).q(2:nzp) - state.fluid(1).q(1:nz))./grid.dzp;
deta2dz = (state.fluid(2).eta(2:nzp) - state.fluid(2).eta(1:nz))./grid.dzp;
dq2dz = (state.fluid(2).q(2:nzp) - state.fluid(2).q(1:nz))./grid.dzp;

% Buoyancy frequency squared
% Note there are several ways to define this, depending on what is assumed
% about the parcel and the environment through which it moves. In all cases
% assume that the pressure profile remains steady.

% Quantities needed in all versions
drdz1 = eos.drdetap1.*deta1dz + eos.drdqp1.*dq1dz + eos.drdp1.*dpdzbar;
drdz2 = eos.drdetap2.*deta2dz + eos.drdqp2.*dq2dz + eos.drdp2.*dpdzbar;

% Version 1
% Assume fluid i rises while fluid j remains steady
% (Note this doesn't make much sense for fluid i=1 when sigma2 is small)
eos.nsq1 = -constants.phys.gravity*eos.sigma2.*eos.rho1.*(drdz2 - eos.drdp1.*dpdzbar);
eos.nsq2 = -constants.phys.gravity*eos.sigma1.*eos.rho2.*(drdz1 - eos.drdp2.*dpdzbar);

% Version 2
% Assume fluid i rises while fluid j falls at velocity
% w_j = -sigma_i w_i / sigma_j
eos.nsq1 = -constants.phys.gravity*eos.rho1.*
           (drdz2 - (eos.sigma1.*eos.drdp2 + eos.sigma2.*eos.drdp1).*dpdzbar);
eos.nsq2 = -constants.phys.gravity*eos.rho2.*
           (drdz1 - (eos.sigma2.*eos.drdp1 + eos.sigma1.*eos.drdp2).*dpdzbar);

% Version 3 (as used in find_eos)
% Consider a parcel that conserves its eta and q moving relative to a
% background in which eta1, eta2, q1, q2, and p are steady.
env_term = eos.sigma1.*drdz1 + eos.sigma2.*.drdz2;;
eos.nsq1 = -constants.phys.gravity*eos.rho1.*(env_term - eos.drdp1.*dpdzbar);
eos.nsq2 = -constants.phys.gravity*eos.rho2.*(env_term - eos.drdp2.*dpdzbar);

% Version 4
% Consider only the properties of fluid i (a parcel of fluid i rises
% relative to the rest of fluid i, ignoring effect of fluid j on p profile)
eos.nsq1 = -constants.phys.gravity*eos.rho1.*(eos.drdetap1.*deta1dz + eos.drdqp1.*dq1dz);
eos.nsq2 = -constants.phys.gravity*eos.rho2.*(eos.drdetap2.*deta2dz + eos.drdqp2.*dq2dz);



