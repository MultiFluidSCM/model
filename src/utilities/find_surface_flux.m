function [ surface_flux ] = find_surface_flux( state, grid, eos, force, psurf, dpdzbar, constants )
% Determine surface fluxes of various quantities and their partitioning
% between fluid 1 and fluid 2


% Fractions of fluxes to assign to each fluid
% Mean surface N^2
% nsq_sfc = eos.sigma1(1)*eos.nsq1(1) + eos.sigma2(1)*eos.nsq2(1)
% Normalized sigma2 to ensure fractions sum to exactly 1
sig2norm = eos.sigma2(1)/(eos.sigma1(1) + eos.sigma2(1));
% if false
if force.sshf > 0
% if nsq_sfc < 0
    % ff2 = sig2norm*(4 - 3*sig2norm);
    % ff2 = 4*sig2norm/(1 + 3*sig2norm);
    % ff2 = 0.4;
    ff2 = max(0.4,min(eos.sigma2(1) + 0.1,1));
    ff1_top = max(0.4,min(eos.sigma1(end) + 0.1,1));
else
    ff2 = sig2norm;
    ff1_top = 1-eos.sigma2(end)/(eos.sigma1(end) + eos.sigma2(end));
end
ff1 = 1 - ff2;
ff2_top = 1 - ff1_top;

% ff1 = 0*eos.sigma2(1) + 1;
% ff2 = 0*eos.sigma2(1);
% ff1_top = 0*eos.sigma2(1);
% ff2_top = 0*eos.sigma2(1) + 1;

% Distribute water flux between the two fluids
surface_flux.q1 = force.sqf*ff1;
surface_flux.q2 = force.sqf*ff2;
surface_flux.q1_top = force.tqf*ff1_top;
surface_flux.q2_top = force.tqf*ff2_top;

% And partial fluxes at first interior level
surface_flux.L1q1 = 0*surface_flux.q1*grid.belowp(1);
surface_flux.L1q2 = 0*surface_flux.q2*grid.belowp(1);
surface_flux.L1q1_top = 0*surface_flux.q1_top*grid.belowp(end);
surface_flux.L1q2_top = 0*surface_flux.q2_top*grid.belowp(end);


% Effective temperature for input of entropy
% Version 1: if we include fluxes at first interior level
T1 = state.fluid(1).T(1);
% T2 = state.fluid(1).T(2);
% m1 = state.fluid(1).m(1);
% m2 = state.fluid(1).m(2);
% mbar = grid.abover(2)*m2 + grid.belowr(2)*m1;
% Teff1 = grid.belowp(1)*(grid.abover(2)*m2*T2 + grid.belowr(2)*m1*T1)/mbar ...
%       + grid.abovep(1)* grid.abover(1)*T1;
% Version 2: if we do not include fluxes at first interior level
Teff1 = T1;
Teff1_top = state.fluid(1).T(end);
 
a = 1-state.fluid(1).q(1);
[g,gp,gt,ga,gpp,gpt,gpa,gtt,gta,gaa,gtwv] = ...
    gibbsav(psurf,Teff1,a,constants.therm);
a_top = 1-state.fluid(1).q(end);
[g_top,gp_top,gt_top,ga_top,gpp_top,gpt_top,gpa_top,gtt_top,gta_top,gaa_top,gtwv_top] = ...
    gibbsav(psurf,Teff1,a_top,constants.therm);

% Surface entropy flux (fluid 1)
surface_flux.eta1 = (force.sshf/Teff1)*ff1 - gtwv*surface_flux.q1;
surface_flux.eta1_top = (force.tshf/Teff1_top)*ff1_top - gtwv_top*surface_flux.q1_top;

% And partial flux at first interior level
surface_flux.L1eta1     = 0*surface_flux.eta1*grid.belowp(1);
surface_flux.L1eta1_top = 0*surface_flux.eta1_top*grid.belowp(end);

% Surface energy flux (for diagnostics)
surface_flux.E1     = Teff1*surface_flux.eta1 + (g - a*ga)*surface_flux.q1;
surface_flux.E1_top = Teff1_top*surface_flux.eta1_top + (g_top - a_top*ga_top)*surface_flux.q1_top;

% Effective temperature for input of entropy
% Version 1: if we include fluxes at first interior level
T1 = state.fluid(2).T(1);
% T2 = state.fluid(2).T(2);
% m1 = state.fluid(2).m(1);
% m2 = state.fluid(2).m(2);
% mbar = grid.abover(2)*m2 + grid.belowr(2)*m1;
% Teff2 = grid.belowp(1)*(grid.abover(2)*m2*T2 + grid.belowr(2)*m1*T1)/mbar ...
%      + grid.abovep(1)* grid.abover(1)*T1;
% Version 2: if we do not include fluxes at first interior level 
Teff2 = T1;
Teff2_top = state.fluid(2).T(end);

a = 1-state.fluid(2).q(1);
[g,gp,gt,ga,gpp,gpt,gpa,gtt,gta,gaa,gtwv] = ...
    gibbsav(psurf,Teff2,a,constants.therm);
a_top = 1-state.fluid(2).q(end);
[g_top,gp_top,gt_top,ga_top,gpp_top,gpt_top,gpa_top,gtt_top,gta_top,gaa_top,gtwv_top] = ...
    gibbsav(psurf,Teff2,a_top,constants.therm);

asat = findasatl(psurf, Teff2, 1-a, constants.therm);
RH = relhum(psurf,Teff2,1-a,constants.therm);
if RH >= 1
    disp(' !! Saturated surface air !!')
    %pause
end

% Surface entropy flux (fluid 2)
surface_flux.eta2     = (force.sshf/Teff2)*ff2 - gtwv*surface_flux.q2;
surface_flux.eta2_top = (force.tshf/Teff2_top)*ff2_top - gtwv_top*surface_flux.q2_top;

% And partial flux at first interior level
surface_flux.L1eta2 = 0*surface_flux.eta2*grid.belowp(1);
surface_flux.L1eta2_top = 0*surface_flux.eta2_top*grid.belowp(end);

% Surface energy flux (for diagnostics)
surface_flux.E2 = Teff2*surface_flux.eta2 + (g - a*ga)*surface_flux.q2;
surface_flux.E2_top = Teff2_top*surface_flux.eta2_top + (g_top - a_top*ga_top)*surface_flux.q2_top;

% ------

% Surface sensible heat flux (bottom and top)
surface_flux.sensible     = force.sshf;
surface_flux.sensible_top = force.tshf;

% Total surface entropy flux
surface_flux.eta = surface_flux.eta1 + surface_flux.eta2;
surface_flux.eta_top = surface_flux.eta1_top + surface_flux.eta2_top;

% Total surface water flux
surface_flux.q = surface_flux.q1 + surface_flux.q2;
surface_flux.q_top = surface_flux.q1_top + surface_flux.q2_top;

% Total surface energy flux
surface_flux.E = surface_flux.E1 + surface_flux.E2;
surface_flux.E_top = surface_flux.E1_top + surface_flux.E2_top;

% Surface buoyancy flux (m^2/s^3)
surface_flux.buoyancy = -constants.phys.gravity*((eos.drdeta1(1)*surface_flux.eta1 ...
                                                + eos.drdq1(1)  *surface_flux.q1) ...
                                               + (eos.drdeta2(1)*surface_flux.eta2 ...
                                                + eos.drdq2(1)  *surface_flux.q2));
surface_flux.buoyancy_top = -constants.phys.gravity*((eos.drdeta1(end)*surface_flux.eta1_top ...
                                                    + eos.drdq1(end)  *surface_flux.q1_top) ...
                                                   + (eos.drdeta2(end)*surface_flux.eta2_top ...
                                                    + eos.drdq2(end)  *surface_flux.q2_top));

end

