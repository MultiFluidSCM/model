function sigma2 = initial_sigma( z_in, initial )
%INITIAL_THETA Set initial potential temperature
%   Provides a potential temperature as a function of height z

n = numel(initial.z);

i = 2;
while z_in > initial.z(i) && i < n
    i = i + 1;
end
sigma2 = initial.sigma2(i-1) + (initial.sigma2(i) - initial.sigma2(i-1)) .* (z_in - initial.z(i-1))./(initial.z(i) - initial.z(i-1));

end