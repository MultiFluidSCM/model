function [Fq, dFqdqa, dFqdqb, dqudz, qubar] = tracer_flux_p( grid, F, q )

% For tracers on p levels:
% compute tracer flux given tracer mixing ratio q
% and mass flux F. Also compute terms needed for
% linearization

nz = grid.nz;
nzp = nz+1;
dzp = grid.dzp;
abovep = grid.abovep;
belowp = grid.belowp;

% For testing, use first order upwind
% Select upwind q's
% Make sure surface flux is zero even if
% F(1) > 0 (due to surface moisture flux)
ixu = [1:nzp] - (F > 0);
ixu(1) = 1;
ixu(nzp) = nz;
qu = q(ixu);
qu(1) = 0;

% Construct fluxes
Fq = F.*qu;

% Derivative of flux wrt q's above and below
dFqdqa = min(F,0);
dFqdqb = max(F,0);

% Average of upwind q's (note reversed weighting)
qubar(1:nz) = abovep(1:nz).*qu(1:nz) + belowp(1:nz).*qu(2:nzp);

% Vertical derivative of upwind q's
dqudz(1:nz) = (qu(2:nzp) - qu(1:nz))./dzp(1:nz);
