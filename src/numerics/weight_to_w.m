function mbar = weight_to_w( grid, m )

% Remap a field from p to w levels using `reversed'
% weighting for conservation

nz = grid.nz;
nzp = nz + 1;

mbar(1:nz) = grid.abover(1:nz).*m;
mbar(nzp) = 0;
mbar(2:nzp) = mbar(2:nzp) + grid.belowr(2:nzp).*m;


end

