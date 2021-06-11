% Preliminary calculations

% Compatibility corrections to ensure older scripts can still run the model
compatibility

if not(settings.restart_simulation)
    % Either set up model for an initial run
    initial_run
else
    % or do a bit-reproducible restart
    read_restart
end

% ---

% Set up timing information
time = settings.time;
settings.time.tstart = current_time;

% First guess for next time step is the current state
state_new = state_old;

% Compute some global diagnostics
gdiags = global_diags(grid,state_new,constants);

% Compute tendencies to obtain diagnostic quantities for plotting
old_diff.flag = 0;
[tend,relabel,eos,force,scales,surface_flux,budgets,work] = ...
         tendencies(grid,state_new,settings,time.t,time.dt,switches,old_diff);

% and set some fields that otherwise would not be set
budgets.w2.wfix   = zeros(1,grid.nz + 1);
budgets.eta1.wfix = zeros(1,grid.nz + 1);
budgets.eta2.wfix = zeros(1,grid.nz + 1);
budgets.q1.wfix   = zeros(1,grid.nz + 1);
budgets.q2.wfix   = zeros(1,grid.nz + 1);

% Output initial state
plottype = 0;
plot_diagnostics

