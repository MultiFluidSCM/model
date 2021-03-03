function [I0, I1, I2] = I_uniform(xs)
%Partial moments of uniform pdf


xstar = sqrt(3);
r2xstar = 0.5/xstar;

if xs < - xstar
    I0 = 1;
    I1 = 0;
    I2 = 1;
elseif xs > xstar
    I0 = 0;
    I1 = 0;
    I2 = 0;
else
    I0 = (xstar - xs)*r2xstar;
    I1 = 0.5*(xstar^2 - xs^2)*r2xstar;
    I2 = (xstar^3 - xs^3)*r2xstar/3;
end


end

