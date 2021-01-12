% Set up model for an initial run


% Set up computational grid
grid = setup_grid();

% Set physical and parameterization constants
constants = set_constants(grid);

% Decide which approximations to impose
switches = set_approximations();

% Initialize accumulated forcing and budget diagnostics
[accum_force,accum] = ini_accum_force();

% Set initial data
state_old = set_initial(grid,constants);

% Set initial time
current_time = 0;

