% Set up model for an initial run


% Set up computational grid
grid = setup_grid();

% Get constants from test case settings
constants = settings.constants;
if ~exist('constants.param.mix')
    constants.param.sort = constants.param;
    constants.param.dwdz = constants.param;
    constants.param.mix = constants.param;
    constants.param.instab = constants.param;
end

% Decide which approximations to impose
switches = set_approximations();

% Initialize accumulated forcing and budget diagnostics
[accum_force,accum] = ini_accum_force();

% Set initial data
state_old = set_initial(grid,constants);

% Set initial time
current_time = 0;

