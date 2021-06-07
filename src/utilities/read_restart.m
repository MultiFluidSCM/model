% Read data needed to do a bit-reproducible restart

% Read data needed to specify the model setup, current state,
% and time series diagnostic data
load(settings.restart_simulation_name);

% Update time settings so the simulation starts and ends for the correct time range.
settings.time.t = current_time;
settings.time.tstop = current_time + settings.time.dt*settings.time.nstop;