% Build linear system for eta variances

% Detrainment factor to allow for effect of variance on detrained value
detfac = M12bar.*relabel.f_sort_chi_hat./max(0.001,sqrt(state_new.fluid(2).vareta));

% Factors needed to allow for buoyancy correlation terms in
% linearization
deta1dz = (eta1(2:nzp) - eta1(1:nz))./grid.dzp;
deta2dz = (eta2(2:nzp) - eta2(1:nz))./grid.dzp;
temp1 = constants.phys.gravity*t_scale1.*deta1dz.*eos.sigma1;
temp2 = constants.phys.gravity*t_scale2.*deta2dz.*eos.sigma2;
temp1_bar = weight_to_w(grid,temp1);
temp2_bar = weight_to_w(grid,temp2);

% Initialize
cc = zeros(5,2*nzp);

ix  = 1:2:2*nzp - 1;
ixb = 3:2:2*nzp - 1;
ixt = 1:2:2*nzp - 3;
% Relabelling term
cc(3 ,ix) = cc(3 ,ix) + M12bar;
cc(4 ,ix) = cc(4 ,ix) - M12bar - (relabel.etahat12 - eta1).*detfac;
% Dissipation term
cc(3 ,ix) = cc(3 ,ix) + m1bar./T_turb1_bar;
% Buoyancy correlation term
%tempap = abovep.*temp1.*eos.rho_deriv_eta1(2:nzp);
%tempbp = belowp.*temp1.*eos.rho_deriv_eta1(1:nz );
%cc(1 ,ixb) = cc(1 ,ixb) - belowr(2:nzp).*tempbp;
%cc(3 ,ixb) = cc(3 ,ixb) - belowr(2:nzp).*tempap;
%cc(3 ,ixt) = cc(3 ,ixt) - abover(1:nz ).*tempbp;
%cc(5 ,ixt) = cc(5 ,ixt) - abover(1:nz ).*tempap;
%cc(3 ,ix) = cc(3 ,ix) - temp1_bar.*eos.rho_deriv_eta1;

ix  = 2:2:2*nzp;
ixb = 4:2:2*nzp;
ixt = 2:2:2*nzp - 2;
% Relabelling term
cc(3 ,ix) = cc(3 ,ix) + M21bar + (relabel.etahat12 - eta2).*detfac;
cc(2 ,ix) = cc(2 ,ix) - M21bar;
% Dissipation term
cc(3 ,ix) = cc(3 ,ix) + m2bar./T_turb2_bar;
% Buoyancy correlation term
%tempap = abovep.*temp2.*eos.rho_deriv_eta2(2:nzp);
%tempbp = belowp.*temp2.*eos.rho_deriv_eta2(1:nz );
%cc(1 ,ixb) = cc(1 ,ixb) - belowr(2:nzp).*tempbp;
%cc(3 ,ixb) = cc(3 ,ixb) - belowr(2:nzp).*tempap;
%cc(3 ,ixt) = cc(3 ,ixt) - abover(1:nz ).*tempbp;
%cc(5 ,ixt) = cc(5 ,ixt) - abover(1:nz ).*tempap;
%cc(3 ,ix) = cc(3 ,ix) - temp2_bar.*eos.rho_deriv_eta2;
