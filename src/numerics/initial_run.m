% Set up model for an initial run


% Set up computational grid
grid = settings.grid;

% Get constants from test case settings
constants = settings.constants;
if not(isfield(constants.param, 'mix'))
    disp("Individual settings for transfer terms not detected");
    constants.param.sort = constants.param;
    constants.param.dwdz = constants.param;
    constants.param.mix = constants.param;
    constants.param.instab = constants.param;
end
disp([num2str(constants.param.mix.bentrainw), " ", num2str(constants.param.mix.bdetrainw)]);

% Decide which approximations to impose
switches = settings.switches;

% Initialize accumulated forcing and budget diagnostics
[accum_force,accum] = ini_accum_force();

% Set initial data
state_old = set_initial(grid, settings);

% Set initial time
current_time = 0;

