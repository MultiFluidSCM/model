function eta = thetae2eta( thetae, q, therm )
%THETAE2ETA Convert equivalent potential temperature to entropy
% Given the wet equivalent potential temperature thetae and
% mass fraction of total water q, determine the entropy eta


% Gibbs functions for liquid water and dry air
[g,gp,gt,ga,gpp,gpt,gpa,gtt,gta,gaa,gtw] = gibbsav(therm.p0d,thetae,1,therm);
[gl,glp,glt,glpp,glpt,gltt] = gibbsl(therm.p0d,theta,therm);
eta = -(1 - q)*gt - q*glt;


end