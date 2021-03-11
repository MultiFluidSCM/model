function [ L_turb ] = find_lturb( grid, nsq, tke, tke_min )
% Compute a turbulence length scale

power = 1;
% power = 3;

nz = numel(grid.zp);

% First length scale is distance from the ground
RL1 = 1.0./grid.zp;

% Second is based on distance a parcel can travel against stratification
RL2 = sqrt(max(nsq,0)./max(tke,tke_min));

% Bound L2 so that dL/dz does not exceed 1 in magnitude
for k = nz-1:-1:1
    RL2(k) = max(RL2(k),RL2(k+1)/(1 + RL2(k+1)*grid.dzw(k+1)));
end
for k = 2:1:nz
    RL2(k) = max(RL2(k),RL2(k-1)/(1 + RL2(k-1)*grid.dzw(k)));
end

% Blend
L_turb = (RL1.^power + RL2.^power).^(-1/power);


% ** Note that, because of the bounding of L2, the dependence of
% L_turb on tke is non-local and is therefore not linearized.



end

