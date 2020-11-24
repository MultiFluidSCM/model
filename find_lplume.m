function [ L_plume ] = find_lplume( grid, nsq, tke, tke_min )
% Compute a plume length scale
% Similar to find_lturb, but do not include distance from the ground
% as the plume is large scale and quasi-laminar near the ground

power = 1;

nz = numel(grid.zp);

% Based on distance a parcel can travel against stratification
RL2 = sqrt(max(nsq,0)./max(tke,tke_min));

% Bound L2 so that dL/dz does not exceed 1 in magnitude
for k = nz-1:-1:1
    RL2(k) = max(RL2(k),RL2(k+1)/(1 + RL2(k+1)*grid.dzw(k+1)));
end
for k = 2:1:nz
    RL2(k) = max(RL2(k),RL2(k-1)/(1 + RL2(k-1)*grid.dzw(k)));
end

% Final expression
L_plume = 1./RL2;


% ** Note that, because of the bounding of L2, the dependence of
% L_plume on tke is non-local and is therefore not linearized.



end

