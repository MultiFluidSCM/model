function [G] = G_uniform(xs)
%Uniform pdf

xstar = sqrt(3);
r2xstar = 0.5/xstar;

if xs < - xstar
    G = 0;
elseif xs > xstar
    G = 0;
else
    G = r2xstar;
end



end


