function force = set_forcing(t)

% Set information needed for forcing

% Set time (x) and sensible (shf) - latent heat fluxes (lhf) in W m-2

% These are for testing/debugging
% x   = [0;  7200; 14400; 27000; 36000; 45000; 52200];
% shf = [-30; -30;     0;     0;   100;   -10; -10];
% lhf = [30;   30;   450;   500;   420;   180;   0];

% These are the correct ARM values
x   = [0;   14400; 23400; 27000; 36000; 45000; 52200];
shf = [-30; 90;    140;   140;   100;   -10;   -10];
lhf = [5; 250;   450;   500;   420;   180;   0];


% Surface fluxes (W / m^2)
for i=1:length(x)-1
   if (x(i) <= t && t < x(i+1))
       force.sshf =   shf(i) + (shf(i+1) - shf(i)).*(t-x(i))./(x(i+1) - x(i)) ;
       force.sqf  = ( lhf(i) + (lhf(i+1) - lhf(i)).*(t-x(i))./(x(i+1) - x(i)) )./2.5e6;
   end
end
if t >=x(end)
   force.sshf = shf(end);
   force.sqf = lhf(end);
end

% force.sshf = 0.0
% force.sqf = 0.0

% Surface moisture flux (kg / m^2 / s)
% Multiply by 2.5e6 to get a rough `latent heat flux'
% But note the concept of `latent heat flux' involves
% several approximations
% force.sqf = 0*10.0e-5;

% Geostrophic wind
force.ug = 10.0;
force.vg = 0.0;


end

