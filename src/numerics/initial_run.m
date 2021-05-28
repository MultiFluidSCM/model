% Set up model for an initial run


% Set up computational grid
grid = settings.grid;

% Get constants from test case settings
constants = settings.constants;

% Decide which approximations to impose
switches = settings.switches;

% Initialize accumulated forcing and budget diagnostics
[accum_force,accum] = ini_accum_force();

% Set initial data
state_old = set_initial(grid, settings);

% Set initial time
current_time = 0;

% Compatibility corrections to ensure older scripts can still run the model
compatibility