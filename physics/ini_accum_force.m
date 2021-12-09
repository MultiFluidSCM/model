function [ accum_force, accum ] = ini_accum_force( )
%INI_ACCUM_FORCE Initialize diagnostics of accumulated forcing
% and energy budgets

% Surface fluxes
accum_force.smf   = 0;  % Mass flux
accum_force.sEf   = 0;  % Energy flux
accum_force.sqf   = 0;  % Moisture flux
accum_force.setaf = 0;  % Entropy flux

% Resolved KE budgets
accum.ke1_pg      = 0;;
accum.ke1_coriol  = 0;
accum.ke1_diff    = 0;
accum.ke1_drag    = 0;
accum.ke1_relabel = 0;
accum.ke2_pg      = 0;;
accum.ke2_coriol  = 0;
accum.ke2_diff    = 0;
accum.ke2_drag    = 0;
accum.ke2_relabel = 0;

% TKE budgets
accum.tke1_shear   = 0;
accum.tke1_buoy    = 0;
accum.tke1_drag    = 0;
accum.tke1_diss    = 0;
accum.tke1_relabel = 0;
accum.tke1_fix     = 0;
accum.tke2_shear   = 0;
accum.tke2_buoy    = 0;
accum.tke2_drag    = 0;
accum.tke2_diss    = 0;
accum.tke2_relabel = 0;
accum.tke2_fix     = 0;

% Entropy sources due to buoyancy flux and turbulent dissipation
accum.eta_bflux    = 0;
accum.eta_dissn    = 0;

end

