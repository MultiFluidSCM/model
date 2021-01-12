function [G] = G_gauss(xs)
%Gaussian pdf

rr2pi = 1/sqrt(2*pi);
G = rr2pi*exp(-xs*xs/2);

end

