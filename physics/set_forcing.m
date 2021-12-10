function force = set_forcing(grid, forcing, t)

% Forcings which aren't time dependent

% Geostrophic wind and radiative cooling
for k = 1:length(grid.zp)
    force.ug(k)   = initial_field(grid.zp(k), forcing.ug.z, forcing.ug.f);
    force.vg(k)   = initial_field(grid.zp(k), forcing.vg.z, forcing.vg.f);
end

% Subsidence
for k = 1:length(grid.zw)
    force.wsub(k) = initial_field(grid.zw(k), forcing.wsub.z, forcing.wsub.f);
    force.rad(k)  = initial_field(grid.zw(k), forcing.rad.z,  forcing.rad.f);
    force.q(k)    = initial_field(grid.zw(k), forcing.q.z,    forcing.q.f);
end

% Time-dependent forcings
force.sshf = interpolate_forcing(forcing.sshf, t);
force.sqf  = interpolate_forcing(forcing.slhf, t)./2.5e6;
force.tshf = interpolate_forcing(forcing.tshf, t);
force.tqf  = interpolate_forcing(forcing.tlhf, t)./2.5e6;

% Note on surface moisture flux (kg / m^2 / s):
% Multiply by 2.5e6 to get a rough `latent heat flux'
% But note the concept of `latent heat flux' involves
% several approximations

end

