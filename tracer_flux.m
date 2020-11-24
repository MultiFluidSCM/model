function [Fq, dFqdqa, dFqdqb, dqudz, qubar] = tracer_flux( grid, Fbar, q )

% For tracers on w levels:
% compute tracer flux given tracer mixing ratio q
% and mass flux Fbar. Also compute terms needed for
% linearization

nz = grid.nz;
nzp = nz+1;
dzw = grid.dzw;
abovew = grid.abovew;
beloww = grid.beloww;

% For testing, use first order upwind
% Select upwind q's
ixu = [1:nz] + (Fbar < 0);
qu = q(ixu);

% Construct fluxes
Fq = Fbar.*qu;

% Derivative of flux wrt q's above and below
dFqdqa = min(Fbar,0);
dFqdqb = max(Fbar,0);

% Average of upwind q's (note reversed weighting)
qubar(1)    = abovew(1)*    q(1)       + beloww(1)*    qu(1);
qubar(2:nz) = abovew(2:nz).*qu(1:nz-1) + beloww(2:nz).*qu(2:nz);
qubar(nzp)  = abovew(nzp)*  qu(nz)     + beloww(nzp)*  q(nzp);

% Vertical derivative of upwind q's
dqudz(1)    = (qu(1)    - q(1))/dzw(1);
dqudz(2:nz) = (qu(2:nz) - qu(1:nz-1))./dzw(2:nz);
dqudz(nzp)  = (q(nzp)   - qu(nz))/dzw(nzp);
