function thetae = eta2thetae( eta, q, T, therm )
%ETA2THETAE Convert entropy to equivalent potential temperature
%   Given entropy eta and total water mass fraction q,
%   compute the equivalent potential temperature
%   by condensing all the water and moving to reference
%   pressure p0d while conserving entropy.
%   T is input as a plausible first guess.


% Newton iteration
tt = T;
for iter = 1:4
    [g,gp,gt,ga,gpp,gpt,gpa,gtt,gta,gaa,gtw] = gibbsav(therm.p0d,tt,1,therm);
    [gl,glp,glt,glpp,glpt,gltt] = gibbsl(therm.p0d,tt,therm);
    res = eta + (1.0d0 - q)*gt + q*glt;
    t_inc = -res/((1.0d0 - q)*gtt + q*gltt);
    tt = tt + t_inc;
end

thetae = tt;



end

