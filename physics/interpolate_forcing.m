% Forcings which are time dependent
function force = interpolate_forcing(forcing, t)

if t <= forcing.t(1)
    force = forcing.f(1);
elseif t >= forcing.t(end)
    force = forcing.f(end);
else
    % Interpolate forcing and fluxes from values provided
    for i=1:length(forcing.t)-1
        if (forcing.t(i) <= t && t < forcing.t(i+1))
            % Surface fluxes (W / m^2)
            force =   forcing.f(i) + (forcing.f(i+1) - forcing.f(i)) ...
                .*(t-forcing.t(i))./(forcing.t(i+1) - forcing.t(i));
        end
    end
end

end