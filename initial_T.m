function T = initial_T( z )
%INITIAL_T Set initial temperature
%   Provides temperature as a function of height z

gamma = -7e-3;
T = 297.2 + gamma*z;

end