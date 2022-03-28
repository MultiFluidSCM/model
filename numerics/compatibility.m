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

% June 2021: Added option for transfers to use pdf properties instead of b-coefficients
if not(isfield(settings.constants.param, 'use_pdf'))
    settings.constants.param.use_pdf = false;
end

% March 2022: New settings for the number of quasi-Newton iterations at the start of the run
if not(isfield(settings, 'solver'))
    settings.solver = struct();
    if not(isfield(settings.solver, 'qn_long_timesteps'))
        settings.solver.qn_start_timesteps  = 0;
        % Number of iterations for thefirst few timesteps (based on qn_long_timesteps)
        settings.solver.qn_iter_max_start   = 8;
        % Number of iterations for future timesteps, once things have stabilised
        settings.solver.qn_iter_max_default = 4;
    end
end