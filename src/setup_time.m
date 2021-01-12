function time = setup_time(current_time)

% Set up timing information, including semi-implicit
% off-centring

% Time to stop
% 3 hours = 10800;
% 6 hours = 21600;
% 7 hours = 25200;
% 9 hours = 32400;
% 11 hours = 39600;
% 12 hours = 43200;
% 14.5 hours = 52200;
time.tstop = 52200;

% Time step
time.dt = 30.0;

% Number of steps to take
time.nstop = round((time.tstop - current_time)/time.dt);
%time.nstop = 1;

% Current time
time.t = current_time;


% Off-centring parameters
time.alpha = 0.55;
time.beta = 1 - time.alpha;


