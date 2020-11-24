function th = initial_theta( z_in )
%INITIAL_THETA Set initial potential temperature
%   Provides a potential temperature as a function of height z

% ARM initial theta profile
% zn      = [0;     50;    350;   650;    700;   1300;   2500;  4400];
% theta_n = [299.0; 301.5; 302.5; 303.53; 303.7; 307.13; 314.0; 330.2];
zn      = [0;     50;    350;   650;    700;   1300;   2500;  5500];
theta_n = [299.0; 301.5; 302.5; 303.53; 303.7; 307.13; 314.0; 343.2];
%theta_n = [300.0; 300.0; 300.0; 300.; 300.0; 300.; 300.0; 300.0];
%disp('*** isentropic ***')
n = numel(zn);

% Unstable CBL profile used in marginal stability paper
%load th_z_dom.mat
%n = numel(z);
%dz = double(z(2) - z(1));
%zn = double(z) - 0.5*dz;
%theta_n = double(thz); % - thetas;


% Interpolate to height z
% for i=1:length(zn)-1
%    if (zn(i) <= z_in && z_in < zn(i+1))
%         th = theta_n(i) + (theta_n(i+1) - theta_n(i)) .* (z_in - zn(i))./(zn(i+1) - zn(i)) ;
%    end
% end
% if z_in >= zn(end)
%    th = theta_n(end);
% end

i = 2;
while z_in > zn(i) && i < n
    i = i + 1;
end
th = theta_n(i-1) + (theta_n(i) - theta_n(i-1)) .* (z_in - zn(i-1))./(zn(i) - zn(i-1));



% Simple profile with constant theta gradient
% Theta gradient
%gamma0 = 1.95e-3;
%gamma0 = 2.93e-3;
%gamma0 = 3.90e-3;
%th = 297.2 + gamma0*z_in;
% th = max(th,297.6);

end