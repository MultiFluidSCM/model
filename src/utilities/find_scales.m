function [ scales ] = find_scales( grid, state, eos, surface_flux, constants, switches )

% Determine characteristic scales

% First unpack a few fields
nz = grid.nz;
nzp = nz + 1;
p    = state.p;
m1   = state.fluid(1).m;
m2   = state.fluid(2).m;
w2   = state.fluid(2).w;
eta1 = state.fluid(1).eta;
eta2 = state.fluid(2).eta;
q1   = state.fluid(1).q;
q2   = state.fluid(2).q;
u1   = state.fluid(1).u;
u2   = state.fluid(2).u;
v1   = state.fluid(1).v;
v2   = state.fluid(2).v;
tke1 = state.fluid(1).tke;
Tw11 = state.fluid(1).Tw(1);
Tw21 = state.fluid(2).Tw(1);

% Remap mass to w levels
m1bar = weight_to_w(grid,m1);
m2bar = weight_to_w(grid,m2);

% ------

% Diagnose characteristic scales

% Surface density, eta, and q
rho00 = m1bar(1) + m2bar(1);
eta00 = (m1bar(1)*eta1(1) + m2bar(1)*eta2(1))/rho00;
q00   = (m1bar(1)*q1(1)   + m2bar(1)*q2(1))/rho00;
T00   = (m1bar(1)*Tw11    + m2bar(1)*Tw21)/rho00;

% Bottom level wind speed
ulev1 = (m1(1)*u1(1) + m2(1)*u2(1))/(m1(1) + m2(1));
vlev1 = (m1(1)*v1(1) + m2(1)*v2(1))/(m1(1) + m2(1));
speed = sqrt(ulev1*ulev1 + vlev1*vlev1);

% Friction velocity
scales.ustar = constants.phys.k*speed/log((grid.zp(1) + constants.param.zrough)/constants.param.zrough);

% Surface buoyancy flux (m^2/s^3)
Bstar = surface_flux.buoyancy;

% Level of neutral buoyancy of a surface parcel, ignoring condensation
% scales.LNBgas = LNBgas(eta2(1),q2(1),Tw21,grid.zw,eos.rhow1,eos.pbar,constants.therm);
scales.LNBgas = LNBgas(eta00,q00,T00,grid.zw,eos.rhow1,eos.pbar,constants.therm);

% Solve iteratively for zstar and wstar
% wstar = 0;
% for iter = 1:3
%     tke_thresh = 0.1*scales.ustar*scales.ustar + 0.01*wstar*wstar;
%     %tke_thresh = max(0.1*max(tke1),2e-6);
%     z_tketop = TKEtop(grid.zp,tke1,tke_thresh);
%     % Convective vertical velocity scale
%     if Bstar > 0
%         wstar = (z_tketop*Bstar)^(1/3);
%     else
%         wstar = 0;
%     end
% end
% scales.tke_thresh = tke_thresh;

% Depth scale zstar
%option_zstar = 4;
%scales.zstar = findzstar(grid.zw,w2,state.fluid(1).Tw,eos.theta_rho1,eos.theta_rho2,option_zstar);
scales.zstar = max(constants.param.zstar_min,scales.LNBgas);
% scales.zstar = z_tketop;

% Convective vertical velocity scale
if Bstar > 0
    wstar = (scales.zstar*Bstar)^(1/3);
else
    wstar = 0;
end
scales.wstar = wstar;

% Friction velocity
scales.ustar = constants.phys.k*speed/log((grid.zp(1) + constants.param.zrough)/constants.param.zrough);
if wstar > 0
    ustar_by_wstar = scales.ustar/wstar;
else
    ustar_by_wstar = 0;
end
scales.ustar_by_wstar = ustar_by_wstar;


% ** Approximate ** specific potential temperature flux
Hstar = surface_flux.sensible/(rho00*constants.therm.Cpd);

% Specific humidity flux
Qstar = surface_flux.q/rho00;

% Specific entropy flux
ETAstar = surface_flux.eta/rho00;

% Approximate theta scale thetastar (may be used to set detrainment rates)
% and approximate humidity scale and entropy scale for near-surface
% relabelling
if wstar > 0
    scales.thetastar = Hstar/wstar;
    scales.qstar     = Qstar/wstar;
    scales.etastar   = ETAstar/wstar;
else
    scales.thetastar = 0;
    scales.qstar     = 0;
    scales.etastar   = 0;
end


% Estimate turbulence length scales
scales.L_turb1 = find_lturb(grid,eos.nsq1,state.fluid(1).tke,constants.param.tke_min);% .* m1 ./ (m1+m2);
scales.L_turb2 = find_lturb(grid,eos.nsq2,state.fluid(2).tke,constants.param.tke_min);% .* m2 ./ (m1+m2);

% and time scales
scales.T_turb1 = scales.L_turb1./sqrt(state.fluid(1).tke);
scales.T_turb2 = scales.L_turb2./sqrt(state.fluid(2).tke);

% and plume length scales
scales.L_plume = find_lplume(grid,eos.nsq2,state.fluid(2).tke,constants.param.tke_min);

end

