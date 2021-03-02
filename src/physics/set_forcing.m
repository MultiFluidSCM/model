function force = set_forcing(forcing, t)

% Interpolate forcing and fluxes from values provided

% Geostrophic wind
force.ug = forcing.ug;
force.vg = forcing.vg;

% Surface fluxes (W / m^2)
for i=1:length(forcing.t)-1
    if (forcing.t(i) <= t && t < forcing.t(i+1))
        force.sshf =   forcing.shf(i) + (forcing.shf(i+1) - forcing.shf(i)) ...
            .*(t-forcing.t(i))./(forcing.t(i+1) - forcing.t(i)) ;
        force.sqf  = ( forcing.lhf(i) + (forcing.lhf(i+1) - forcing.lhf(i)) ...
            .*(t-forcing.t(i))./(forcing.t(i+1) - forcing.t(i)) )./2.5e6;
    end
end
% Set forcing to last available value if time beyond interpolation range
if t >=forcing.t(end)
    force.sshf = forcing.shf(end);
    force.sqf = forcing.lhf(end);
end

% Note on surface moisture flux (kg / m^2 / s):
% Multiply by 2.5e6 to get a rough `latent heat flux'
% But note the concept of `latent heat flux' involves
% several approximations

end

