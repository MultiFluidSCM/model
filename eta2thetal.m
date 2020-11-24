function thetal = eta2thetal( eta, q, T, therm, p00 )
%ETA2THETAL Convert entropy to liquid water potential temperature
%   Given entropy eta and total water mass fraction q,
%   compute the liquid water potential temperature
%   T is input as a plausible first guess.


% Newton iteration
tt = T;
for iter = 1:4
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a ] = gibbs(p00,tt,q,therm);
    res = eta + gt;
    t_inc = -res/gtt;
    tt = tt + t_inc;
end

thetal = tt;

end

