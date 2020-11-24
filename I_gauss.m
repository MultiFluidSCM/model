function [I0, I1, I2] = I_gauss(xs)
%Partial moments of Gaussian


rr2 = 1/sqrt(2);
rr2pi = 1/sqrt(2*pi);

I0 = 0.5*(1 + erf(-xs*rr2));
I1 = rr2pi*exp(-xs*xs/2);
I2 = xs*I1 + I0;


end

