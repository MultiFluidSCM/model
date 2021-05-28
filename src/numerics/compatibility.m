% Compatibility corrections to ensure older scripts can still run the model

% In January 2021, transfer parameters for each individual entrainment/detrainment type were added
if not(isfield(constants.param, 'mix'))
    disp("Individual settings for transfer terms not detected");
    constants.param.sort = constants.param;
    constants.param.dwdz = constants.param;
    constants.param.mix = constants.param;
    constants.param.instab = constants.param;
end

% In May 2021, extra transfer parameters for mixing in the cloud layer were added
if not(isfield(constants.param, 'mix_cloud'))
    disp("Settings for mixing transfer in cloud not detected");
    constants.param.mix_cloud = constants.param.mix;
end