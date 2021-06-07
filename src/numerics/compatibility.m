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
if not(isfield(settings.constants.param.mix, 'tke_factor'))
    settings.constants.param.mix.tke_factor = 1;
end

% June 2021: Settings to restart/continue a previous simulation are now set in test_cases folders
if not(isfield(settings, 'restart_simulation'))
    settings.restart_simulation = false;
    settings.restart_simulation_name = 'restart_00001740';
end