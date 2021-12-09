function [ rh ] = relhum( p, T, ww, therm )
%RELHUM Find relative humidity p/psat

% First find saturation value of a, the dry air mass fraction in the
% gaseous part.
asat = findasatl(p, T, ww, therm);

if (ww >= 1 - asat)
    % Saturated case
    a = asat;
else
    % Unsaturated case, all water is vapour
    a = (1-ww);
end

% Relative humidity
e = therm.epsilon;
rh = ((1 - a   )*(1 + asat*(e - 1))) ...
   / ((1 - asat)*(1 + a   *(e - 1)));

end

