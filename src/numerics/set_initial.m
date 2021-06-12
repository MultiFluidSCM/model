function state = set_initial( grid, settings )

% Define initial state

% Copy grid arrays for clarity
nz = grid.nz;
nzp = nz + 1;
zw = grid.zw;
zp = grid.zp;
dzw = grid.dzw;
dzp = grid.dzp;
abovew = grid.abovew;
beloww = grid.beloww;
abovep = grid.abovep;
belowp = grid.belowp;
aboves = grid.aboves;
belows = grid.belows;

% Useful constants
constants = settings.constants;
% Gravity
gravity = constants.phys.gravity;

% Surface pressure
%p0s = 100000;
p0s =  97000;


% Specify profiles in balance with a given thetal(z)
disp('** note that initial_theta is used as thetal **')
% pause
hydrostatic_thetal

% or a given T(z)
%hydrostatic_T


% Now loop again and back out rho-level density and temperature
for k = 1:nz
    p00 = p(k);
    ww     = aboves(k)*water(k+1) + belows(k)*water(k);
    etabar = aboves(k)*eta(k+1)   + belows(k)*eta(k);
    tt = abovep(k)*Tw(k+1) + belowp(k)*Tw(k);
    for iter1 = 1:10
        [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(p00,tt,ww,constants.therm);
        r1 = gt + etabar;
        t_inc = - r1/gtt;
        tt = tt + t_inc;
    end
    % Save rho-level density and temperature
    rho(k) = 1/gp;
    T(k) = tt;
end


% Set winds to geostrophic values, which are set in set_forcing
force = set_forcing(settings.forcing, 0);

% Save in state structure
%state.fluid(1).m   = (1.0 - constants.param.confrac)*rho;
state.fluid(1).eta = eta;
state.fluid(1).q   = water;
state.fluid(1).T   = T;
state.fluid(1).Tw  = Tw;
state.fluid(1).w(1:nzp) = 0;
state.fluid(1).u(1:nz) = force.ug;
state.fluid(1).v(1:nz) = force.vg;
state.fluid(1).tke = ones(1,nz)*constants.param.tke_min;
state.fluid(1).vareta(1:nz) = 0;
state.fluid(1).varq(1:nz)   = 0;
state.p            = p;
% Fluid 2 is the same as fluid 1 except for m
state.fluid(2)     = state.fluid(1);
% Initial mass fractions
sigma00 = constants.param.sigma00;
state.fluid(2).m   = sigma00.*rho;
state.fluid(1).m   = rho - state.fluid(2).m;


end

