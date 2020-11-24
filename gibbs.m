function [ g, gp, gt, gw, gpp, gpt, gtt, gpw, gtw, gww, a ] = gibbs( p, T, ww, therm )
%GIBBS Gibbs fn and derivatives for air - water vapour - liquid water mix
%   Given pressure p, temperature T and mass fraction of total water ww,
%   compute the Gibbs thermodynamic potential along with its first and
%   second derivatives.

% First find saturation value of a, the dry air mass fraction in the
% gaseous part.
asat = findasatl(p, T, ww, therm);

if (ww >= 1 - asat)
    % Saturated case
    a = asat;
    wl = (ww - (1-a))/a;
    wg = 1 - wl;
    % Contribution from gaseous part
    [g,gp,gt,ga,gpp,gpt,gpa,gtt,gta,gaa,gtw] = gibbsav(p,T,a,therm);
    % Contribution from liquid part
    [gl,glp,glt,glpp,glpt,gltt] = gibbsl(p,T,therm);
    
    gw = (gl - g)/a;
    LLp = gp - a*gpa - glp;
    LLt = gt - a*gta - glt;
    rden = 1/(a*a*gaa);
    gpw = (glp - gp)/a;
    gtw = (glt - gt)/a;
    g = wg*g + wl*gl;
    gp = wg*gp + wl*glp;
    gt = wg*gt + wl*glt;
    gpp = wg*gpp + wl*glpp - wg*LLp*LLp*rden;
    gpt = wg*gpt + wl*glpt - wg*LLp*LLt*rden;
    gtt = wg*gtt + wl*gltt - wg*LLt*LLt*rden;
    gww = 0;
    
else
    % Unsaturated case, all water is vapour
    a = (1-ww);
    [g,gp,gt,ga,gpp,gpt,gpa,gtt,gta,gaa,gtw] = gibbsav(p,T,a,therm);
    gw = -ga;
    gpw = -gpa;
    gtw = -gta;
    gww = gaa;
end



end

