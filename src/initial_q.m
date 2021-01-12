function qv = initial_q( z )
%INITIAL_Q Set initial water vapour mixing ratio
%   Provides a water vapour as a function of height z

% ARM initial q profile
zn  = [0;       50;       350;      650;     700;     1300;    2500;   5500];
q_n = [15.2e-3; 15.17e-3; 14.98e-3; 14.8e-3; 14.7e-3; 13.5e-3; 3.0e-3; 3.0e-3];
% q_n = [1.e-4;   1.e-4;    1.e-4;    1.e-4;   1.e-4;   1.e-4;   1.e-4;  1.e-4];
% disp('*** dry run ***')
% Convert from mixing ratio to specific humidity
q_n = q_n./(1 + q_n);

for i=1:length(zn)-1
   if (zn(i) <= z && z < zn(i+1))
        qv = q_n(i)+(q_n(i+1)-q_n(i)) .* (z-zn(i))./(zn(i+1)-zn(i)) ;
   end
end
if z >= zn(end)
    qv = q_n(end);
end


% Simple profile with exponential decay
% Surface specific humidity
% ww0 = 0.0160
% Scale height for specific humidity for a simple
% exponentially decaying profile
% hww = 2000.0;
%qv = ww0*exp(-z/hww);

% Dry case
% qv = 0;

end