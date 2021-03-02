% Plot a range of diagnostics

% Font size
fs = 16;

% Units for z axis on most plots
% km
zunitsp = grid.zpkm;
zunitsw = grid.zwkm;
zplottop = 4.4;
labelz = 'z (km)';
% Or non-dimensionalized
% zunitsp = grid.zp/scales.zstar;
% zunitsw = grid.zw/scales.zstar;
% zplottop = 1.5;
% labelz = 'z / z_*';


% Unpack fields
nz = grid.nz;
nzp = nz + 1;
dzp = grid.dzp;
dzw = grid.dzw;
abovep = grid.abovep;
belowp = grid.belowp;
abovew = grid.abovew;
beloww = grid.beloww;
abover = grid.abover;
belowr = grid.belowr;
m1   = state_new.fluid(1).m;
m2   = state_new.fluid(2).m;
sigma1 = m1./eos.rho1;
sigma2 = m2./eos.rho2;
w1   = state_new.fluid(1).w;
w2   = state_new.fluid(2).w;
eta1 = state_new.fluid(1).eta;
eta2 = state_new.fluid(2).eta;
q1   = state_new.fluid(1).q;
q2   = state_new.fluid(2).q;
p    = state_new.p;
T1   = state_new.fluid(1).T;
T2   = state_new.fluid(2).T;
Tw1  = state_new.fluid(1).Tw;
Tw2  = state_new.fluid(2).Tw;
u1   = state_new.fluid(1).u;
u2   = state_new.fluid(2).u;
v1   = state_new.fluid(1).v;
v2   = state_new.fluid(2).v;
tke1 = state_new.fluid(1).tke;
tke2 = state_new.fluid(2).tke;
m1bar = work.m1bar;
m2bar = work.m2bar;
F1bar = work.F1bar;
F2bar = work.F2bar;
Feta1 = work.Feta1;
Feta2 = work.Feta2;
Deta1 = work.Deta1;
Deta2 = work.Deta2;
Fq1 = work.Fq1;
Fq2 = work.Fq2;
Dq1 = work.Dq1;
Dq2 = work.Dq2;
mbar = m1bar + m2bar;
zstar = scales.zstar;
wstar = scales.wstar;
etastar = scales.etastar;
thetastar = scales.thetastar;
force2 = set_forcing(settings.forcing, time.t);
qstar = scales.qstar;
setaf = surface_flux.eta;
sqf   = force2.sqf;

% Average to p levels
eta1bar = abovep.*eta1(2:nzp)...
        + belowp.*eta1(1:nz );
q1bar   = abovep.*q1(2:nzp)...
        + belowp.*q1(1:nz );
w1bar   = abovep.*w1(2:nzp)...
        + belowp.*w1(1:nz );
eta2bar = abovep.*eta2(2:nzp)...
        + belowp.*eta2(1:nz );
q2bar   = abovep.*q2(2:nzp)...
        + belowp.*q2(1:nz );
w2bar   = abovep.*w2(2:nzp)...
        + belowp.*w2(1:nz );
    
% Average to w levels
sigma1w = weight_to_w(grid,sigma1);
sigma2w = weight_to_w(grid,sigma2);

% ----------

% Diagnose cloud fractions
diagnose_cloud_frac

% Append current time to array for plotting time series
if exist('ts')
    plotstep = numel(ts.time) + 1;
else
    plotstep = 1;
end
ts.time(plotstep) = time.t;

if settings.switches.plot

    % Figure 1: Basic fields mass frac, w, eta, q, u and v, buoyancy, N^2,
    % thetal, RH, tke
    if settings.switches.plot_basic_fields
        plot_basic_fields
    end
    
    % Figure 2: Profile snapshots and time series diagnostics
    if settings.switches.plot_time_series
        plot_time_series
    end

    % Figure 3: Budget diagnostics
    if settings.switches.plot_budgets
        plot_budgets
    end

    % Figure 4: Turbulence related diagnostics
    if settings.switches.plot_turbulence
        plot_turbulence
    end

    % Figure 18: For development: actual and candidate entrainment and
    % detrainment profiles
    if settings.switches.plot_ED
        plot_ED
    end

    % Figure 19: Plot buoyancy of a selection of parcels lifted adiabatically
    if settings.switches.plot_adiabats
        plot_adiabats
    end

    % Figure 20: Plot estimated standard deviations of updraft tracers
    if settings.switches.plot_std
        plot_std
    end

    % Figure 21: Plot estimated cloud fractions using APDF
    if settings.switches.plot_cloud_frac
        plot_cloud_frac
    end

    % Figure 22: Compare N and turbulence inverse time scale
    if settings.switches.plot_T_turb
        plot_T_turb
    end

    % Figure 23: Plot variance tendencies for updraft tracers
    if settings.switches.plot_var_tend
        plot_var_tend
    end

    % Figure 24: Plot diffusive buoyancy fluxes
    if settings.switches.plot_bflux
        plot_bflux
    end

    % Figure 25: Plot diffusive and buoyancy correlation eta fluxes
    if settings.switches.plot_bc_flux
        plot_bc_flux
    end

    % Figure 26: Check that eta and q variance tendencies balance
    if settings.switches.plot_var_budgets
        plot_var_budgets
    end

    % Figure 27: Plot time series of max residuals
    if settings.switches.plot_res_time_series
        plot_res_time_series
    end

    % Figure 28: 
    if settings.switches.plot_ql_var
        plot_ql_var
    end

end


pause(0.01)