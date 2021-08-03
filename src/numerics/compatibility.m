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
