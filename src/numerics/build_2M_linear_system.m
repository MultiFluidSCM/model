% Build the matrix for second moments system:
% tke - tracer fluxes - variances


% Initialize
dd = zeros(21,10*nz);


% TKE1 equation
ix =  1:10:10*nz-9;
ikb = 2:nzp;
ikt = 1:nz;
ikc = 2:nz;
% Tendency term
dd(11,ix) = dd(11,ix) + m1;
% Transport terms    
dd( 1,ix) = dd( 1,ix) - adt*work.dFtke1dtkeb(ikt)./dzp;
dd(11,ix) = dd(11,ix) - adt*(work.dFtke1dtkea(ikt) - work.dFtke1dtkeb(ikb))./dzp;
dd(21,ix) = dd(21,ix) + adt*work.dFtke1dtkea(ikb)./dzp;
% Diffusion terms
dd( 1,ix) = dd( 1,ix) - adt*work.dDtke1dtkeb(ikt)./dzp;
dd(11,ix) = dd(11,ix) - adt*(work.dDtke1dtkea(ikt) - work.dDtke1dtkeb(ikb))./dzp;
dd(21,ix) = dd(21,ix) + adt*work.dDtke1dtkea(ikb)./dzp;
% Dissipation term
dd(11,ix) = dd(11,ix) + adt*m1.*(1.5./scales.L_turb1).*sqrt(tke1);
% Relabelling terms - only keep dependence on tke1
dd(11,ix) = dd(11,ix) + adt*M21;

% TKE2 equation
ix =  2:10:10*nz-8;
ikb = 2:nzp;
ikt = 1:nz;
ikc = 2:nz;
% Tendency term
dd(11,ix) = dd(11,ix) + m2;
% Transport terms    
dd( 1,ix) = dd( 1,ix) - adt*work.dFtke2dtkeb(ikt)./dzp;
dd(11,ix) = dd(11,ix) - adt*(work.dFtke2dtkea(ikt) - work.dFtke2dtkeb(ikb))./dzp;
dd(21,ix) = dd(21,ix) + adt*work.dFtke2dtkea(ikb)./dzp;
% Diffusion terms
dd( 1,ix) = dd( 1,ix) - adt*work.dDtke2dtkeb(ikt)./dzp;
dd(11,ix) = dd(11,ix) - adt*(work.dDtke2dtkea(ikt) - work.dDtke2dtkeb(ikb))./dzp;
dd(21,ix) = dd(21,ix) + adt*work.dDtke2dtkea(ikb)./dzp;
% Dissipation term
dd(11,ix) = dd(11,ix) + adt*m2.*(1.5./scales.L_turb2).*sqrt(tke2);
% Relabelling terms - only keep dependence on tke2
dd(11,ix) = dd(11,ix) + adt*M12;

% (w-eta)1 correlation
% Dummy equation for the moment
ix =  3:10:10*nz-7;
ikb = 2:nzp;
ikt = 1:nz;
ikc = 2:nz;
% Tendency term
dd(11,ix) = dd(11,ix) + 1.0;

% (w-eta)2 correlation
% Dummy equation for the moment
ix =  4:10:10*nz-6;
ikb = 2:nzp;
ikt = 1:nz;
ikc = 2:nz;
% Tendency term
dd(11,ix) = dd(11,ix) + 1.0;

% (w-q)1 correlation
% Dummy equation for the moment
ix =  5:10:10*nz-5;
ikb = 2:nzp;
ikt = 1:nz;
ikc = 2:nz;
% Tendency term
dd(11,ix) = dd(11,ix) + 1.0;

% (w-q)2 correlation
% Dummy equation for the moment
ix =  6:10:10*nz-4;
ikb = 2:nzp;
ikt = 1:nz;
ikc = 2:nz;
% Tendency term
dd(11,ix) = dd(11,ix) + 1.0;


% Some quantities needed for variance equations

 % Vertical pressure gradient
dpdz(2:nz)   = (p(2:nz) - p(1:nz-1))./dzw(2:nz);
dpdz(1) = dpdz(2);
dpdz(nzp) = dpdz(nz);
dpdzbar = abovep.*dpdz(2:nzp) + belowp.*dpdz(1:nz);
dpdz(1) = 0;
dpdz(nzp) = 0;

% Timescale for dissipation of flux in buoyancy correlation terms
t_scale1 = 1.5*scales.L_turb1./sqrt(tke1);
t_scale2 = 1.5*scales.L_turb2./sqrt(tke2);

% Factors needed to allow for buoyancy correlation terms in
% linearization
% Assume zero correlation between eta and q
deta1dz = max(0,(eta1(2:nzp) - eta1(1:nz))./grid.dzp);
deta2dz = max(0,(eta2(2:nzp) - eta2(1:nz))./grid.dzp);
factor1 = 2*t_scale1.*dpdzbar.*m1.*eos.drdetap1.*deta1dz;
factor2 = 2*t_scale2.*dpdzbar.*m2.*eos.drdetap2.*deta2dz;
    
% eta variance1 equation
ix =  7:10:10*nz-3;
ik = 1:nz;
% Buoyancy correlation term
dd(11,ix) = dd(11,ix) + settings.buoy_correl_eta*factor1(ik);
% Dissipation term
dd(11,ix) = dd(11,ix) + m1(ik)./scales.T_turb1(ik);
% Relabelling terms
dd(11,ix) = dd(11,ix) + M12(ik);
dd(12,ix) = dd(12,ix) - M12(ik);

% eta variance2 equation
ix =  8:10:10*nz-2;
ik = 1:nz;
% Buoyancy correlation term
dd(11,ix) = dd(11,ix) + settings.buoy_correl_eta*factor2(ik);
% Dissipation term
dd(11,ix) = dd(11,ix) + m2(ik)./scales.T_turb2(ik);
% Relabelling terms
dd(11,ix) = dd(11,ix) + M21(ik);
dd(10,ix) = dd(10,ix) - M21(ik);

% *** detfac = 0*M12bar(k)*relabel.f_sort_chi_hat(k)/max(0.001,sqrt(state_new.fluid(2).vareta(k)));
% A11 =  M12(k) + m1(k)/scales.T_turb1(k) + settings.buoy_correl_eta*factor1(k);
% A12 = -M12(k); % *** - (relabel.etahat12(k) - eta1(k))*detfac;
% A21 = -M21(k);
% A22 =  M21(k) + m2(k)/scales.T_turb2(k) + settings.buoy_correl_eta*factor2(k); % *** + (relabel.etahat12(k) - eta2(k))*detfac;

% Factors needed to allow for buoyancy correlation terms in
% linearization
% Assume zero correlation between eta and q
dq1dz = (q1(2:nzp) - q1(1:nz))./grid.dzp;
dq2dz = (q2(2:nzp) - q2(1:nz))./grid.dzp;
factor1 = 2*t_scale1.*dpdzbar.*m1.*eos.drdqp1.*dq1dz;
factor2 = 2*t_scale2.*dpdzbar.*m2.*eos.drdqp2.*dq2dz;

% q variance1 equation
ix =  9:10:10*nz-1;
ik = 1:nz;
% Buoyancy correlation term
dd(11,ix) = dd(11,ix) + settings.buoy_correl_q*factor1(ik);
% Dissipation term
dd(11,ix) = dd(11,ix) + m1(ik)./scales.T_turb1(ik);
% Relabelling terms
dd(11,ix) = dd(11,ix) + M12(ik);
dd(12,ix) = dd(12,ix) - M12(ik);

% q variance2 equation
ix =  10:10:10*nz;
ik = 1:nz;
% Buoyancy correlation term
dd(11,ix) = dd(11,ix) + settings.buoy_correl_q*factor2(ik);
% Dissipation term
dd(11,ix) = dd(11,ix) + m2(ik)./scales.T_turb2(ik);
% Relabelling terms
dd(11,ix) = dd(11,ix) + M21(ik);
dd(10,ix) = dd(10,ix) - M21(ik);

% % detfac = M12bar(k)*relabel.f_sort_chi_hat(k)/max(1e-6,sqrt(state_new.fluid(2).varq(k)));
% A11 =  M12(k) + m1(k)/scales.T_turb1(k) + settings.buoy_correl_q*factor1(k);
% A12 = -M12(k); % *** - (relabel.qhat12(k) - q1(k))*detfac;
% A21 = -M21(k);
% A22 =  M21(k) + m2(k)/scales.T_turb2(k) + settings.buoy_correl_q*factor2(k); % *** + (relabel.qhat12(k) - q2(k))*detfac;






