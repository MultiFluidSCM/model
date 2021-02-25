% Accumulate various budget related quantities

% accdt should be set to bdt when calling from advance and to adt
% when calling from quasinewton

% Unpack some fields
nz = grid.nz;
nzp = nz + 1;
dzp = grid.dzp;
dzw = grid.dzw;
u1   = state_new.fluid(1).u;
u2   = state_new.fluid(2).u;
v1   = state_new.fluid(1).v;
v2   = state_new.fluid(2).v;
w1   = state_new.fluid(1).w;
w2   = state_new.fluid(2).w;

% Surface forcing terms in mass, energy, water, and entropy budgets
accum_force.smf   = accum_force.smf    + accdt*force.sqf;
accum_force.sEf   = accum_force.sEf    + accdt*surface_flux.E; %force.sEf;
accum_force.sqf   = accum_force.sqf    + accdt*force.sqf;
accum_force.setaf = accum_force.setaf  + accdt*surface_flux.eta; %force.setaf;

% Terms in budgets of resolved KE
ke1 = 0.5*(u1.^2 + v1.^2 + grid.aboves.*w1(2:nzp).^2 + grid.belows.*w1(1:nz).^2);
ke2 = 0.5*(u2.^2 + v2.^2 + grid.aboves.*w2(2:nzp).^2 + grid.belows.*w2(1:nz).^2);
accum.ke1_pg      = accum.ke1_pg      + accdt*sum(w1.*tend.fluid(1).mw.pgterm.*dzw);
accum.ke1_coriol  = accum.ke1_coriol  + accdt*sum(u1.*tend.fluid(1).mu.coriolis.*dzp ...
                                                + v1.*tend.fluid(1).mv.coriolis.*dzp);
accum.ke1_diff    = accum.ke1_diff    + accdt*sum(w1.*(tend.fluid(1).mw.diffuse ...
                                                     + tend.fluid(1).mw.diffent).*dzw) ...
                                      + accdt*sum(u1.*(tend.fluid(1).mu.diffuse ...
                                                     + tend.fluid(1).mu.diffent).*dzp ...
                                                + v1.*(tend.fluid(1).mv.diffuse ...
                                                     + tend.fluid(1).mv.diffent).*dzp);
accum.ke1_drag    = accum.ke1_drag    + accdt*sum(w1.*tend.fluid(1).mw.drag.*dzw);
accum.ke1_relabel = accum.ke1_relabel + accdt*sum(w1.*tend.fluid(1).mw.relabel.*dzw) ...
                                      + accdt*sum(u1.*tend.fluid(1).mu.relabel.*dzp ...
                                                + v1.*tend.fluid(1).mv.relabel.*dzp ...
                                                - ke1.*tend.fluid(1).m.relabel.*dzp);
accum.ke2_pg      = accum.ke2_pg      + accdt*sum(w2.*tend.fluid(2).mw.pgterm.*dzw);
accum.ke2_coriol  = accum.ke2_coriol  + accdt*sum(u2.*tend.fluid(2).mu.coriolis.*dzp ...
                                                + v2.*tend.fluid(2).mv.coriolis.*dzp);
accum.ke2_diff    = accum.ke2_diff    + accdt*sum(w2.*(tend.fluid(2).mw.diffuse ...
                                                     + tend.fluid(2).mw.diffent).*dzw) ...
                                      + accdt*sum(u2.*(tend.fluid(2).mu.diffuse ...
                                                     + tend.fluid(2).mu.diffent).*dzp ...
                                                + v2.*(tend.fluid(2).mv.diffuse ...
                                                     + tend.fluid(2).mv.diffent).*dzp);
accum.ke2_drag    = accum.ke2_drag    + accdt*sum(w2.*tend.fluid(2).mw.drag.*dzw);
accum.ke2_relabel = accum.ke2_relabel + accdt*sum(w2.*tend.fluid(2).mw.relabel.*dzw) ...
                                      + accdt*sum(u2.*tend.fluid(2).mu.relabel.*dzp ...
                                                + v2.*tend.fluid(2).mv.relabel.*dzp ...
                                                - ke2.*tend.fluid(2).m.relabel.*dzp);

% Terms in budgets of TKE
accum.tke1_shear   = accum.tke1_shear   + accdt*sum(tend.fluid(1).mtke.shear.*dzp);
accum.tke1_buoy    = accum.tke1_buoy    + accdt*sum(tend.fluid(1).mtke.bflux.*dzp);
accum.tke1_drag    = accum.tke1_drag    + accdt*sum(tend.fluid(1).mtke.drag.*dzp);
accum.tke1_diss    = accum.tke1_diss    + accdt*sum(tend.fluid(1).mtke.dissn.*dzp);
accum.tke1_relabel = accum.tke1_relabel + accdt*sum(tend.fluid(1).mtke.relabel.*dzp);
accum.tke2_shear   = accum.tke2_shear   + accdt*sum(tend.fluid(2).mtke.shear.*dzp);
accum.tke2_buoy    = accum.tke2_buoy    + accdt*sum(tend.fluid(2).mtke.bflux.*dzp);
accum.tke2_drag    = accum.tke2_drag    + accdt*sum(tend.fluid(2).mtke.drag.*dzp);
accum.tke2_diss    = accum.tke2_diss    + accdt*sum(tend.fluid(2).mtke.dissn.*dzp);
accum.tke2_relabel = accum.tke2_relabel + accdt*sum(tend.fluid(2).mtke.relabel.*dzp);


% Entropy source due to buoyancy flux and dissipation
accum.eta_bflux = accum.eta_bflux + accdt*sum(tend.fluid(1).meta.bflux.*dzw ...
                                            + tend.fluid(2).meta.bflux.*dzw);
accum.eta_dissn = accum.eta_dissn + accdt*sum(tend.fluid(1).meta.dissn.*dzw ...
                                            + tend.fluid(2).meta.dissn.*dzw);

