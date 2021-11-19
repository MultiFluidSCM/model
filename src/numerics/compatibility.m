% Compatibility corrections to ensure older scripts can still run the model

% Jan 2021: Transfer parameters for each individual entrainment/detrainment added
if not(isfield(settings.constants.param, 'mix'))
    disp("Individual settings for transfer terms not detected");
    settings.constants.param.sort = settings.constants.param;
    settings.constants.param.dwdz = settings.constants.param;
    settings.constants.param.mix = settings.constants.param;
    settings.constants.param.instab = settings.constants.param;
end

% May 2021: Extra transfer parameters for mixing in the cloud layer added
if not(isfield(settings.constants.param, 'mix_cloud'))
    disp("Settings for mixing transfer in cloud not detected");
    settings.constants.param.mix_cloud = settings.constants.param.mix;
end

% May 2021: Tuning coefficient for mixing entrainment/detrainment for scaling.
% May be removed once suitable coefficient has been established
% June 2021: Separate TKE factor introduced for tke1 and tke2
if not(isfield(settings.constants.param.mix, 'tke1_factor'))
    disp("Settings for mixing TKE factor not detected");
    settings.constants.param.mix.tke1_factor = 0.;
    settings.constants.param.mix.tke2_factor = 0.25;
end

% June 2021: Settings to restart/continue a previous simulation are now set in test_cases folders
if not(isfield(settings, 'restart_simulation'))
    disp("Settings for simulation restart not detected, starting model from t=0");
    settings.restart_simulation = false;
    settings.restart_simulation_name = 'restart_00001740';
end

% June 2021: Backward compatibility for buoyancy correlation term switches
if not(isfield(settings, 'buoy_correl_eta'))
    settings.buoy_correl_eta = 0;
end
if not(isfield(settings, 'buoy_correl_q'))
    settings.buoy_correl_q   = 0;
end

% June 2021: Added parameter to control cloud threshold for cloud base and cloud height diagnostic
if not(isfield(settings.constants.param, 'cld_thresh'))
    disp("Cloud threshold not detected");
    settings.constants.param.cld_thresh = 0.01;
end

if not(isfield(settings.constants.param, 'sigma_weighted_tke'))
    disp("Using non-sigma-weighted TKE for turbulent length scales");
    settings.constants.param.sigma_weighted_tke = false;
end

% June 2021: Added option for transfers to use pdf properties instead of b-coefficients
if not(isfield(settings.constants.param, 'use_pdf'))
    settings.constants.param.use_pdf = false;
end
if not(isfield(settings.constants.param.dwdz, 'use_pdf'))
    disp("Using b-coefficients for dw/dz transfer (if dw/dz transfer enabled)");
    settings.constants.param.dwdz.use_pdf = false;
end
if not(isfield(settings.constants.param.instab, 'use_pdf'))
    disp("Using b-coefficients for instability transfer (if instability transfer enabled)");
    settings.constants.param.instab.use_pdf = false;
end
if not(isfield(settings.constants.param.mix, 'use_pdf'))
    disp("Using b-coefficients for mixing transfer (if mixing transfer enabled)");
    settings.constants.param.mix.use_pdf = false;
end

% July 2021: Mellor-Yamada-Nakanishi-Niino constants included in turbulence
% modelling. The values required to reproduce the previous model behaviour
% are as follows.
if not(isfield(settings.constants.param, 'MYNN'))
    settings.constants.param.MYNN.A1 = 1e20;
    settings.constants.param.MYNN.A2 = sqrt(0.5);
    settings.constants.param.MYNN.B1 = sqrt(8);
    settings.constants.param.MYNN.B2 = sqrt(8);
    settings.constants.param.MYNN.C1 = 0;
    settings.constants.param.MYNN.C2 = 0;
    settings.constants.param.MYNN.C3 = 0;
    settings.constants.param.MYNN.C4 = 0;
    settings.constants.param.MYNN.C5 = 0;
end

% August 2021: Added option to choose output times for LES comparison
if not(isfield(settings, 'output_times'))
    settings.output_times = [14000, 21800, 32600, 42800];
end

% September 2021: Removed initial theta and rv profiles so they are no longer hard-coded
if not(isfield(settings, 'initial_theta'))
    disp("Using default initial theta profiles for the ARM case");
    settings.initial_theta.z     = [    0;    50;   350;    650;   700;   1300;  2500;  5500];
    settings.initial_theta.theta = [299.0; 301.5; 302.5; 303.53; 303.7; 307.13; 314.0; 343.2];
end
if not(isfield(settings, 'initial_rv'))
    disp("Using default initial rv profiles for the ARM case");
    settings.initial_rv.z  = [      0;       50;      350;      650;    700;    1300;   2500;   5500];
    settings.initial_rv.rv = [15.2e-3; 15.17e-3; 14.98e-3; 14.8e-3; 14.7e-3; 13.5e-3; 3.0e-3; 3.0e-3];
end
if not(isfield(settings, 'initial_qv'))
    disp("Converting initial rv to initial qv");
    settings.initial_qv.z  = settings.initial_rv.z;
    settings.initial_qv.qv = settings.initial_rv.rv./(1 + settings.initial_rv.rv);
end
if not(isfield(settings, 'initial_sigma'))
    disp("Using default initial sigma profiles for the ARM case");
    settings.initial_sigma.z      = [0; 1];
    settings.initial_sigma.sigma2 = [settings.constants.param.sigma00; settings.constants.param.sigma00];
end

% September 2021: Added option to add fluxes at the top of the domain
if not(isfield(settings.forcing, 'tshf'))
    disp("Using no prescribed fluxes at the top of the domain");
    settings.forcing.tshf = 0*settings.forcing.t;
    settings.forcing.tlhf = 0*settings.forcing.t;
end

% October 2021: New forcing terms used in BOMEX case
if not(isfield(settings.forcing, 'ug_z'))
    settings.forcing.ug_z = 0;
    settings.forcing.vg_z = 0;
end

% October 2021: Added initial proviles for the horizontal velocities
if not(isfield(settings, 'initial_u'))
    disp("Using default initial u profiles for the ARM case");
    settings.initial_u.z = [0];
    settings.initial_u.u = [10.];
end
if not(isfield(settings, 'initial_v'))
    disp("Using default initial v profiles for the ARM case");
    settings.initial_v.z = [0];
    settings.initial_v.v = [0];
end

% October 2021: New subsidence forcing term used in BOMEX case
if not(isfield(settings.forcing, 'wsub'))
    disp("Using default subsidence of 0");
    settings.forcing.wsub_z = 0;
    settings.forcing.wsub   = 0;
end

% October 2021: New moisture forcing
if not(isfield(settings.forcing, 'q'))
    disp("Using default moisture forcing of 0");
    settings.forcing.q_z = 0;
    settings.forcing.q   = 0;
end

% October 2021: New radiative cooling forcing term used in BOMEX case
if not(isfield(settings.forcing, 'rad'))
    disp("Using default radiative cooling of 0");
    settings.forcing.rad_z = 0;
    settings.forcing.rad   = 0;
end

% October 2021: Option to set the surface pressure added
if not(isfield(settings, 'surface_pressure'))
    disp("Using default ARM surface pressure of 97000");
    settings.surface_pressure = 97000;
end

% November 2021: 
if not(isfield(settings.constants.param.instab, 'entrain_factor'))
    settings.constants.param.instab.entrain_factor = 0.2;
end
if not(isfield(settings.constants.param.instab, 'detrain_factor'))
    settings.constants.param.instab.detrain_factor = 0.2;
end