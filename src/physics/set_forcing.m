function force = set_forcing(grid, forcing, t)

% Forcings which aren't time dependent

% Geostrophic wind and radiative cooling
for k = 1:length(grid.zp)
    force.ug(k)   = initial_field(grid.zp(k), forcing.ug_z, forcing.ug);
    force.vg(k)   = initial_field(grid.zp(k), forcing.vg_z, forcing.vg);
end

% Subsidence
for k = 1:length(grid.zw)
    force.wsub(k) = initial_field(grid.zw(k), forcing.wsub_z, forcing.wsub);
    force.rad(k)  = initial_field(grid.zw(k), forcing.rad_z,  forcing.rad);
    force.q(k)    = initial_field(grid.zw(k), forcing.q_z,    forcing.q);
end

% Forcings which are time dependent
% Interpolate forcing and fluxes from values provided
for i=1:length(forcing.t)-1
    if (forcing.t(i) <= t && t < forcing.t(i+1))
        % Surface fluxes (W / m^2)
        force.sshf =   forcing.shf(i) + (forcing.shf(i+1) - forcing.shf(i)) ...
            .*(t-forcing.t(i))./(forcing.t(i+1) - forcing.t(i)) ;
        force.sqf  = ( forcing.lhf(i) + (forcing.lhf(i+1) - forcing.lhf(i)) ...
            .*(t-forcing.t(i))./(forcing.t(i+1) - forcing.t(i)) )./2.5e6;
        
        % Fluxes at the top of the domain (W / m^2)
        force.tshf =   forcing.tshf(i) + (forcing.tshf(i+1) - forcing.tshf(i)) ...
            .*(t-forcing.t(i))./(forcing.t(i+1) - forcing.t(i)) ;
        force.tqf  = ( forcing.tlhf(i) + (forcing.tlhf(i+1) - forcing.tlhf(i)) ...
            .*(t-forcing.t(i))./(forcing.t(i+1) - forcing.t(i)) )./2.5e6;
    end
end

% Set forcing to last available value if time beyond interpolation range
if t >=forcing.t(end)
    force.sshf = forcing.shf(end);
    force.sqf  = forcing.lhf(end)./2.5e6;
    force.tshf = forcing.tshf(end);
    force.tqf  = forcing.tlhf(end)./2.5e6;
end

% Note on surface moisture flux (kg / m^2 / s):
% Multiply by 2.5e6 to get a rough `latent heat flux'
% But note the concept of `latent heat flux' involves
% several approximations

end

