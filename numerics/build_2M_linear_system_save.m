% Build the matrix for second moments system:
% tke - tracer fluxes - variances

% Useful preliminary quantities

 % Vertical pressure gradient
dpdz(2:nz)   = (p(2:nz) - p(1:nz-1))./dzw(2:nz);
dpdz(1) = dpdz(2);
dpdz(nzp) = dpdz(nz);
dpdzbar = aboves.*dpdz(2:nzp) + belows.*dpdz(1:nz);
dpdz(1) = 0;
dpdz(nzp) = 0;

% Vertical derivatives of eta and q
deta1dz = (eta1(2:nzp) - eta1(1:nz))./grid.dzp;
deta2dz = (eta2(2:nzp) - eta2(1:nz))./grid.dzp;
dq1dz   = (q1(2:nzp)   - q1(1:nz)  )./grid.dzp;
dq2dz   = (q2(2:nzp)   - q2(1:nz)  )./grid.dzp;

% Timescale for dissipation of flux in buoyancy correlation terms
% Now in work.T_sflux
% t_scale1 = 1.5*scales.L_turb1./sqrt(tke1);
% t_scale2 = 1.5*scales.L_turb2./sqrt(tke2);

% Factors appearing in buoyancy correlation terms
factor1_eta = work.T_sflux1.*dpdzbar.*m1.*eos.drdetap1;
factor1_q   = work.T_sflux1.*dpdzbar.*m1.*eos.drdqp1;
factor2_eta = work.T_sflux2.*dpdzbar.*m2.*eos.drdetap2;
factor2_q   = work.T_sflux2.*dpdzbar.*m2.*eos.drdqp2;

% Derivative of diffusivity wrt tke at p levels.
% Assumes all diffusivities are the same
dKdtke1 = (scales.dLdtke1 + 0.5*scales.L_turb1./tke1).*sqrt(tke1);
dKdtke2 = (scales.dLdtke2 + 0.5*scales.L_turb2./tke2).*sqrt(tke2);

% Derivative of flux timescale wrt tke at p levels.
dTdtke1 = 3*sqrt(0.5)*settings.constants.param.MYNN.A2*(scales.dLdtke1 - 0.5*scales.L_turb1./tke1)./sqrt(tke1);
dTdtke2 = 3*sqrt(0.5)*settings.constants.param.MYNN.A2*(scales.dLdtke2 - 0.5*scales.L_turb2./tke2)./sqrt(tke2);

% ----

% Initialize
dd = zeros(25,12*nz);


% TKE1 equation
ix =  1:12:12*nz-11;
ikb = 2:nzp;
ikt = 1:nz;
ikc = 2:nz;
% Tendency term
dd(13,ix) = dd(13,ix) + m1;
% Transport terms    
dd( 1,ix) = dd( 1,ix) - adt*work.dFtke1dtkeb(ikt)./dzp;
dd(13,ix) = dd(13,ix) - adt*(work.dFtke1dtkea(ikt) - work.dFtke1dtkeb(ikb))./dzp;
dd(25,ix) = dd(25,ix) + adt*work.dFtke1dtkea(ikb)./dzp;
% Diffusion terms
dd( 1,ix) = dd( 1,ix) - adt*work.dDtke1dtkeb(ikt)./dzp;
dd(13,ix) = dd(13,ix) - adt*(work.dDtke1dtkea(ikt) - work.dDtke1dtkeb(ikb))./dzp;
dd(25,ix) = dd(25,ix) + adt*work.dDtke1dtkea(ikb)./dzp;
% Buoyancy flux term
dd(15,ix) = dd(15,ix) - adt*dpdzbar.*eos.drdetap1;
dd(17,ix) = dd(17,ix) - adt*dpdzbar.*eos.drdqp1;
% Dissipation term
dd(13,ix) = dd(13,ix) + adt*m1.*work.dissn_rate_tke1.*work.rate_lin_fac1;
% Relabelling terms - only keep dependence on tke1
dd(13,ix) = dd(13,ix) + adt*M21;

% TKE2 equation
ix =  2:12:12*nz-10;
ikb = 2:nzp;
ikt = 1:nz;
ikc = 2:nz;
% Tendency term
dd(13,ix) = dd(13,ix) + m2;
% Transport terms    
dd( 1,ix) = dd( 1,ix) - adt*work.dFtke2dtkeb(ikt)./dzp;
dd(13,ix) = dd(13,ix) - adt*(work.dFtke2dtkea(ikt) - work.dFtke2dtkeb(ikb))./dzp;
dd(25,ix) = dd(25,ix) + adt*work.dFtke2dtkea(ikb)./dzp;
% Diffusion terms
dd( 1,ix) = dd( 1,ix) - adt*work.dDtke2dtkeb(ikt)./dzp;
dd(13,ix) = dd(13,ix) - adt*(work.dDtke2dtkea(ikt) - work.dDtke2dtkeb(ikb))./dzp;
dd(25,ix) = dd(25,ix) + adt*work.dDtke2dtkea(ikb)./dzp;
% Buoyancy flux term
dd(15,ix) = dd(15,ix) - adt*dpdzbar.*eos.drdetap2;
dd(17,ix) = dd(17,ix) - adt*dpdzbar.*eos.drdqp2;
% Dissipation term
dd(13,ix) = dd(13,ix) + adt*m2.*work.dissn_rate_tke2.*work.rate_lin_fac2;
% Relabelling terms - only keep dependence on tke2
dd(13,ix) = dd(13,ix) + adt*M12;


% Fluid 1 SG eta flux
ix =  3:12:12*nz-9;
% Diagonal term
dd(13,ix) = dd(13,ix) + 1.0;
% Downgradient diffusion
dd(11,ix) = dd(11,ix) + m1.*deta1dz.*dKdtke1;
% Buoyancy correlation term
if settings.buoy_correl_eta
    dd(17,ix) = dd(17,ix) - factor1_eta;
    dd(11,ix) = dd(11,ix) - m1.*dpdzbar.*eos.drdetap1.*state_new.fluid(1).vareta.*dTdtke1;
    dd(21,ix) = dd(21,ix) - factor1_q;
end

% Fluid 2 SG eta flux
ix =  4:12:12*nz-8;
% Diagonal term
dd(13,ix) = dd(13,ix) + 1.0;
% Downgradient diffusion
dd(11,ix) = dd(11,ix) + m2.*deta2dz.*dKdtke2;
% Buoyancy correlation term
if settings.buoy_correl_eta
    dd(17,ix) = dd(17,ix) - factor2_eta;
    dd(11,ix) = dd(11,ix) - m2.*dpdzbar.*eos.drdetap2.*state_new.fluid(2).vareta.*dTdtke2;
    dd(21,ix) = dd(21,ix) - factor2_q;
end

% Fluid 1 SG q flux
ix =  5:12:12*nz-7;
% Diagonal term
dd(13,ix) = dd(13,ix) + 1.0;
% Downgradient diffusion
dd( 9,ix) = dd( 9,ix) + m1.*dq1dz.*dKdtke1;
% Buoyancy correlation term
if settings.buoy_correl_q
    dd(17,ix) = dd(17,ix) - factor1_q;
    dd( 9,ix) = dd( 9,ix) - m1.*dpdzbar.*eos.drdqp1.*state_new.fluid(1).varq.*dTdtke1;
    dd(19,ix) = dd(19,ix) - factor1_eta;
end

% Fluid 2 SG q flux
ix =  6:12:12*nz-6;
% Diagonal term
dd(13,ix) = dd(13,ix) + 1.0;
% Downgradient diffusion
dd( 9,ix) = dd( 9,ix) + m2.*dq2dz.*dKdtke2;
% Buoyancy correlation term
if settings.buoy_correl_q
    dd(17,ix) = dd(17,ix) - factor2_q;
    dd( 9,ix) = dd( 9,ix) - m2.*dpdzbar.*eos.drdqp2.*state_new.fluid(2).varq.*dTdtke2;
    dd(19,ix) = dd(19,ix) - factor2_eta;
end

    
% eta variance1 equation
ix =  7:12:12*nz-5;
ik = 1:nz;
if settings.buoy_correl_eta
    % Buoyancy correlation term
    dd(13,ix) = dd(13,ix) + 2*factor1_eta.*work.deta1dz_modified;
    dd(17,ix) = dd(17,ix) + 2*factor1_q  .*work.deta1dz_modified;
% Quasi-diffusion terms to allow for feedback via deta/dz
    % - *** Borrow tke coefficients for testing ***
    dd( 1,ix) = dd( 1,ix) - adt*work.dissn_rate_var1.*work.dDtke1dtkeb(ikt)./dzp;
    dd(13,ix) = dd(13,ix) - adt*work.dissn_rate_var1.*(work.dDtke1dtkea(ikt) - work.dDtke1dtkeb(ikb))./dzp;
    dd(25,ix) = dd(25,ix) + adt*work.dissn_rate_var1.*work.dDtke1dtkea(ikb)./dzp;
end
% Dissipation term
dd(13,ix) = dd(13,ix) + m1(ik).*work.dissn_rate_var1(ik);
% Relabelling terms
dd(13,ix) = dd(13,ix) + M12(ik);
dd(14,ix) = dd(14,ix) - M12(ik);

% eta variance2 equation
ix =  8:12:12*nz-4;
ik = 1:nz;
if settings.buoy_correl_eta
    % Buoyancy correlation term
    dd(13,ix) = dd(13,ix) + 2*factor2_eta.*work.deta2dz_modified;
    dd(17,ix) = dd(17,ix) + 2*factor2_q  .*work.deta2dz_modified;
    % Quasi-diffusion terms to allow for feedback via deta/dz
    % - *** Borrow tke coefficients for testing ***
    dd( 1,ix) = dd( 1,ix) - adt*work.dissn_rate_var2.*work.dDtke2dtkeb(ikt)./dzp;
    dd(13,ix) = dd(13,ix) - adt*work.dissn_rate_var2.*(work.dDtke2dtkea(ikt) - work.dDtke2dtkeb(ikb))./dzp;
    dd(25,ix) = dd(25,ix) + adt*work.dissn_rate_var2.*work.dDtke2dtkea(ikb)./dzp;
end
% Dissipation term
dd(13,ix) = dd(13,ix) + m2(ik).*work.dissn_rate_var2(ik);
% Relabelling terms
dd(13,ix) = dd(13,ix) + M21(ik);
dd(12,ix) = dd(12,ix) - M21(ik);

% Factors needed to account for effect of q2 variance on sorting detrainment
rrootqvar2 = 1./max(1e-6,sqrt(state_new.fluid(2).varq));
detfac_up = M12.*grid.aboves.*relabel.f_sort_chi_hat(2:nzp);
detfac_dn = M12.*grid.belows.*relabel.f_sort_chi_hat(1:nz);

% q variance1 equation
ix =  9:12:12*nz-3;
ik = 1:nz;
if settings.buoy_correl_q
    % Buoyancy correlation term
    dd(13,ix) = dd(13,ix) + 2*factor1_q  .*work.dq1dz_modified;
    dd(15,ix) = dd(15,ix) + 2*factor1_eta.*work.dq1dz_modified;
    % Quasi-diffusion terms to allow for feedback via dq/dz
    % - *** Borrow tke coefficients for testing ***
    dd( 1,ix) = dd( 1,ix) - adt*work.dissn_rate_var1.*work.dDtke1dtkeb(ikt)./dzp;
    dd(13,ix) = dd(13,ix) - adt*work.dissn_rate_var1.*(work.dDtke1dtkea(ikt) - work.dDtke1dtkeb(ikb))./dzp;
    dd(25,ix) = dd(25,ix) + adt*work.dissn_rate_var1.*work.dDtke1dtkea(ikb)./dzp;
end
% Dissipation term
dd(13,ix) = dd(13,ix) + m1(ik).*work.dissn_rate_var1(ik);
% Relabelling terms
dd(13,ix) = dd(13,ix) + M12(ik);
dd(14,ix) = dd(14,ix) - M12(ik);
% ... including effect of q2 variance on sorting
dd(14,ix) = dd(14,ix) - detfac_up.*(relabel.qhat12(2:nzp) - q1(2:nzp)).*grid.beloww(2:nzp).*rrootqvar2 ...
                      - detfac_dn.*(relabel.qhat12(1:nz ) - q1(1:nz )).*grid.abovew(1:nz ).*rrootqvar2;

% q variance2 equation
ix  = 10:12:12*nz-2;
ixb = 22:12:12*nz-2;
ixt = 10:12:12*nz-14;
ik = 1:nz;
if settings.buoy_correl_q
    % Buoyancy correlation term
    dd(13,ix) = dd(13,ix) + 2*factor2_q  .*work.dq2dz_modified;
    dd(15,ix) = dd(15,ix) + 2*factor2_eta.*work.dq2dz_modified;
    % Quasi-diffusion terms to allow for feedback via dq/dz
    % - *** Borrow tke coefficients for testing ***
    dd( 1,ix) = dd( 1,ix) - adt*work.dissn_rate_var2.*work.dDtke2dtkeb(ikt)./dzp;
    dd(13,ix) = dd(13,ix) - adt*work.dissn_rate_var2.*(work.dDtke2dtkea(ikt) - work.dDtke2dtkeb(ikb))./dzp;
    dd(25,ix) = dd(25,ix) + adt*work.dissn_rate_var2.*work.dDtke2dtkea(ikb)./dzp;
end
% Dissipation term
dd(13,ix) = dd(13,ix) + m2(ik).*work.dissn_rate_var2;
% Relabelling terms
dd(13,ix) = dd(13,ix) + M21(ik);
dd(12,ix) = dd(12,ix) - M21(ik);
% ... including effect of q2 variance on sorting
% dd( 1,ixb) = dd( 1,ixb) + detfac_dn(2:nz  ).*(relabel.qhat12(2:nz ) - q2(2:nz )).*grid.beloww(2:nz  ).*rrootqvar2(1:nz-1);
dd(13,ix ) = dd(13,ix ) + detfac_up        .*(relabel.qhat12(2:nzp) - q2(2:nzp)).*grid.beloww(2:nzp ).*rrootqvar2 ...
                        + detfac_dn        .*(relabel.qhat12(1:nz ) - q2(1:nz )).*grid.abovew(1:nz  ).*rrootqvar2;
% dd(25,ixt) = dd(25,ixt) + detfac_up(1:nz-1).*(relabel.qhat12(2:nz ) - q2(2:nz )).*grid.abovew(1:nz-1).*rrootqvar2(2:nz  );

% eta-q covariance1 equation
ix =  11:12:12*nz-1;
ik = 1:nz;
% Buoyancy correlation term
if settings.buoy_correl_eta
    dd(13,ix) = dd(13,ix) + factor1_q  .*work.dq1dz_modified;
    dd( 9,ix) = dd( 9,ix) + factor1_eta.*work.dq1dz_modified;
end
if settings.buoy_correl_q
    dd(13,ix) = dd(13,ix) + factor1_eta.*work.deta1dz_modified;
    dd(11,ix) = dd(11,ix) + factor1_q  .*work.deta1dz_modified;
end
if settings.buoy_correl_eta | settings.buoy_correl_q
    % Quasi-diffusion terms to allow for feedback via dq/dz, detadz
    % - *** Borrow tke coefficients for testing ***
    dd( 1,ix) = dd( 1,ix) - adt*work.dissn_rate_var1.*work.dDtke1dtkeb(ikt)./dzp;
    dd(13,ix) = dd(13,ix) - adt*work.dissn_rate_var1.*(work.dDtke1dtkea(ikt) - work.dDtke1dtkeb(ikb))./dzp;
    dd(25,ix) = dd(25,ix) + adt*work.dissn_rate_var1.*work.dDtke1dtkea(ikb)./dzp;
end
% Dissipation term
dd(13,ix) = dd(13,ix) + m1(ik).*work.dissn_rate_var1;
% Relabelling terms
dd(13,ix) = dd(13,ix) + M12(ik);
dd(14,ix) = dd(14,ix) - M12(ik);

% eta-q covariance2 equation
ix =  12:12:12*nz;
ik = 1:nz;
% Buoyancy correlation term
if settings.buoy_correl_eta
    dd(13,ix) = dd(13,ix) + factor2_q  .*work.dq2dz_modified;
    dd( 9,ix) = dd( 9,ix) + factor2_eta.*work.dq2dz_modified;
end
if settings.buoy_correl_q
    dd(13,ix) = dd(13,ix) + factor2_eta.*work.deta2dz_modified;
    dd(11,ix) = dd(11,ix) + factor2_q  .*work.deta2dz_modified;
end
if settings.buoy_correl_eta | settings.buoy_correl_q
    % Quasi-diffusion terms to allow for feedback via dq/dz, detadz
    % - *** Borrow tke coefficients for testing ***
    dd( 1,ix) = dd( 1,ix) - adt*work.dissn_rate_var2.*work.dDtke2dtkeb(ikt)./dzp;
    dd(13,ix) = dd(13,ix) - adt*work.dissn_rate_var2.*(work.dDtke2dtkea(ikt) - work.dDtke2dtkeb(ikb))./dzp;
    dd(25,ix) = dd(25,ix) + adt*work.dissn_rate_var2.*work.dDtke2dtkea(ikb)./dzp;
end
% Dissipation term
dd(13,ix) = dd(13,ix) + m2(ik).*work.dissn_rate_var2(ik);
% Relabelling terms
dd(13,ix) = dd(13,ix) + M21(ik);
dd(12,ix) = dd(12,ix) - M21(ik);



