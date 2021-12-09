function eta = thetal2eta( thetal, q, therm, p00 )
%THETAL2ETA Convert liquid water potential temperature to entropy
% Given the liquid water potential temperature thetal and
% mass fraction of total water q, determine the entropy eta


% Gibbs functions for liquid water and dry air
[g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a ] = gibbs(p00,thetal,q,therm);
eta = -gt;


end