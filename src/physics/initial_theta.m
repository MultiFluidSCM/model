function th = initial_theta( z_in, initial )
%INITIAL_THETA Set initial potential temperature
%   Provides a potential temperature as a function of height z

n = numel(initial.z);

i = 2;
while z_in > initial.z(i) && i < n
    i = i + 1;
end
th = initial.theta(i-1) + (initial.theta(i) - initial.theta(i-1)) .* (z_in - initial.z(i-1))./(initial.z(i) - initial.z(i-1));

end