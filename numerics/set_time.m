function time = set_time(time_settings)

% Set up timing information, including semi-implicit
% off-centring

% 3 hours = 10800;
% 6 hours = 21600;
% 7 hours = 25200;
% 9 hours = 32400;
% 11 hours = 39600;
% 12 hours = 43200;
% 14.5 hours = 52200;
time.tstart = time_settings.tstart;
time.tstop = time_settings.tstop;

% Time step
time.dt = time_settings.dt;

% Number of steps to take
time.nstop = round((time.tstop - time.tstart)/time.dt);

% Current time
time.t = time.tstart;

% Off-centring parameters
time.alpha = 0.55;
time.beta = 1 - time.alpha;

end