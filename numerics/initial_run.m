% Set up model for an initial run

% Set up the discretisation in time
time = set_time(settings.time);

% Set up computational grid
grid = set_grid(settings.grid);

% Get constants from test case settings
constants = settings.constants;

% Compute profile of geopotential
constants.phi = constants.phys.gravity * grid.zp;

% Decide which approximations to impose
switches = settings.switches;

% Initialize accumulated forcing and budget diagnostics
[accum_force,accum] = ini_accum_force();

% Set initial data
state_old = set_initial(grid, settings);

% Set initial time
current_time = 0;

