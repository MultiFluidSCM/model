function [ L_turb, dLdtke ] = find_lturb_deriv( grid, nsq, tke, tke_min )

% Compute a turbulence length scale
% Also estimate its local derivative wrt tke. This is complicated because
% of the bounding of L2, but empirical tests suggest it captures the
% derivative reasonably well except in some places where L2 is bounded.

power = 1;
% power = 3;

% Avoid divisions by zero
tkex = max(tke,tke_min);

nz = numel(grid.zp);

% First length scale is distance from the ground
RL1 = 1.0./grid.zp;

% Second is based on distance a parcel can travel against stratification
RL20 = sqrt(max(nsq,0)./tkex);

% Bound L2 so that dL/dz does not exceed 1 in magnitude
RL2(nz) = RL20(nz);
for k = nz-1:-1:1
    RL2(k) = max(RL20(k),RL2(k+1)/(1 + RL2(k+1)*grid.dzw(k+1)));
end
for k = 2:1:nz
    RL2(k) = max(RL2(k),RL2(k-1)/(1 + RL2(k-1)*grid.dzw(k)));
end

% Blend
L_turb = (RL1.^power + RL2.^power).^(-1/power);

% Compute derivative
dLdtke = 0.5*(L_turb.*(L_turb.*RL20).^power)./tkex;



end

