% Compatibility corrections to ensure older scripts can still run the model

% In January 2021, transfer parameters for each individual entrainment/detrainment type were added
if not(isfield(settings.constants.param, 'mix'))
    disp("Individual settings for transfer terms not detected");
    settings.constants.param.sort = settings.constants.param;
    settings.constants.param.dwdz = settings.constants.param;
    settings.constants.param.mix = settings.constants.param;
    settings.constants.param.instab = settings.constants.param;
end

% In May 2021, extra transfer parameters for mixing in the cloud layer were added
if not(isfield(settings.constants.param, 'mix_cloud'))
    disp("Settings for mixing transfer in cloud not detected");
    settings.constants.param.mix_cloud = settings.constants.param.mix;
end