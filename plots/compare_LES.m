% Save data for direct comparison with LES results

% List of times at which to compare with LES
SCM_output_times = settings.output_times;
% settings.output_times = [14000:600:42800];
nst = numel(SCM_output_times); 

% Time step size
dt = time.dt;


% Is it time to save data?
lsave = 0;
ist = 0;
while ist < nst & ~lsave
    ist = ist + 1;
    st = SCM_output_times(ist);
    if (time.t - dt < st & time.t >= st)
        lsave = 1;
        SCM_times(ist) = st;
    end
end

% If required then store data
if lsave
    
    % Unpack some fields
    m1   = state_new.fluid(1).m;
    m2   = state_new.fluid(2).m;
    sigma1 = m1./eos.rho1;
    sigma2 = m2./eos.rho2;
    sigma1w = weight_to_w(grid,sigma1);
    sigma2w = weight_to_w(grid,sigma2);
    w1   = state_new.fluid(1).w;
    w2   = state_new.fluid(2).w;
    % eta1 = state_new.fluid(1).eta;
    % eta2 = state_new.fluid(2).eta;
    q1   = state_new.fluid(1).q;
    q2   = state_new.fluid(2).q;
    p    = state_new.p;
    % T1   = state_new.fluid(1).T;
    % T2   = state_new.fluid(2).T;
    Tw1  = state_new.fluid(1).Tw;
    Tw2  = state_new.fluid(2).Tw;
    u1   = state_new.fluid(1).u;
    u2   = state_new.fluid(2).u;
    v1   = state_new.fluid(1).v;
    v2   = state_new.fluid(2).v;
    tke1 = state_new.fluid(1).tke;
    tke2 = state_new.fluid(2).tke;
    m1bar = work.m1bar;
    m2bar = work.m2bar;
    Cpd = constants.therm.Cpd;
    
    % Determine liquid water
    disp('compare_LES: use eos.ql1, eos.ql2 for liquid water')
    for k = 1:nzp
        if k == 1
            pbar = grid.extrapb1*p(1) + grid.extrapb2*p(2);
        elseif k == nzp
            pbar = grid.extraptnz*p(nz) + grid.extraptnzm*p(nz-1);
        else
            pbar   = grid.abovew(k)*p(k) ...
                   + grid.beloww(k)*p(k-1);
        end
        [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(pbar,Tw1(k),q1(k),constants.therm);
        vapour1(k) = (1 - q1(k))*(1 - a)/a;
        liquid1(k) = q1(k) - vapour1(k);
        [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(pbar,Tw2(k),q2(k),constants.therm);
        vapour2(k) = (1 - q2(k))*(1 - a)/a;
        liquid2(k) = q2(k) - vapour2(k);
    end
    
    % TKE
    % KE of mean resolved
    umean = (m1.*u1 + m2.*u2)./(m1 + m2);
    vmean = (m1.*v1 + m2.*v2)./(m1 + m2);
    wmean = (m1bar.*w1 + m2bar.*w2)./(m1bar + m2bar);
    % Resolved w variance
    Ruu1 = (umean - u1).^2;
    Ruu2 = (umean - u2).^2;
    Rvv1 = (vmean - v1).^2;
    Rvv2 = (vmean - v2).^2;
    Rww1 = (wmean - w1).^2;
    Rww2 = (wmean - w2).^2;
    % Sub-filter w variance
    uu1 = (2/3)*state_new.fluid(1).tke;
    uu2 = (2/3)*state_new.fluid(2).tke;
    vv1 = (2/3)*state_new.fluid(1).tke;
    vv2 = (2/3)*state_new.fluid(2).tke;
    ww1 = (2/3)*state_new.fluid(1).tke;
    ww2 = (2/3)*state_new.fluid(2).tke;
    % Total w variance per unit mass
    uutot = (m1.*(Ruu1 + uu1) + m2.*(Ruu2 + uu2))./(m1 + m2);
    vvtot = (m1.*(Rvv1 + vv1) + m2.*(Rvv2 + vv2))./(m1 + m2);
    wwtot = (m1.*(grid.aboves.*Rww1(2:nzp) + grid.belows.*Rww1(1:nz) + ww1) + m2.*(grid.aboves.*Rww2(2:nzp) + grid.belows.*Rww2(1:nz) + ww2))./(m1 + m2);
    % Resolved TKE per unit mass
    Rtke1 = 0.5*((umean - u1).^2 + (vmean - v1).^2 + grid.aboves.*Rww1(2:nzp) + grid.belows.*Rww1(1:nz));
    Rtke2 = 0.5*((umean - u2).^2 + (vmean - v2).^2 + grid.aboves.*Rww2(2:nzp) + grid.belows.*Rww2(1:nz));
    % Total TKE per unit mass
    tketot = (m1.*(Rtke1 + tke1) + m2.*(Rtke2 + tke2))./(m1 + m2);
    
    % Resolved q variance
    qmean = (m1bar.*q1 + m2bar.*q2)./(m1bar + m2bar);
    qvmean = (m1bar.*vapour1 + m2bar.*vapour2)./(m1bar + m2bar);
    qlmean = (m1bar.*liquid1 + m2bar.*liquid2)./(m1bar + m2bar);
    Rqq1 = (q1 - qmean).^2;
    Rqq2 = (q2 - qmean).^2;
    
    % Resolved theta(_l) variance.
    thetamean = (m1bar.*eos.theta1 + m2bar.*eos.theta2)./(m1bar + m2bar);
    Rthth1 = (eos.theta1 - thetamean).^2;
    Rthth2 = (eos.theta2 - thetamean).^2;
    
    % Pack data ready for output
    
    % Mean profiles
    SCM_m(:,ist) = m1 + m2;
    SCM_m_1(:,ist) = m1;
    SCM_m_2(:,ist) = m2;
    SCM_mw(:,ist) = m1bar + m2bar;
    SCM_m_1w(:,ist) = m1bar;
    SCM_m_2w(:,ist) = m2bar;
    SCM_rho(:,ist) = (m1.*eos.rho1 + m2.*eos.rho2)./(m1 + m2);
    SCM_rho_1(:,ist) = eos.rho1;
    SCM_rho_2(:,ist) = eos.rho2;
    SCM_rhow(:,ist) = weight_to_w(grid, (m1.*eos.rho1 + m2.*eos.rho2)./(m1 + m2));
    SCM_rho_1w(:,ist) = weight_to_w(grid, eos.rho1);
    SCM_rho_2w(:,ist) = weight_to_w(grid, eos.rho2);
    
    % Mass fractions at p levels and w levels
    SCM_sigma1(:,ist) = sigma1;
    SCM_sigma2(:,ist) = sigma2;
    SCM_sigma1w(:,ist) = sigma1w;
    SCM_sigma2w(:,ist) = sigma2w;
    
    % Mass flux
    SCM_mf(:,ist) = m1bar.*w1 + m2bar.*w2;
    SCM_mf_1(:,ist) = m1bar.*w1;
    SCM_mf_2(:,ist) = m2bar.*w2;
    
    % Horizontal velocity
    SCM_u(:,ist) = umean;
    SCM_u_1(:,ist) = u1;
    SCM_u_2(:,ist) = u2;
    SCM_v(:,ist) = vmean;
    SCM_v_1(:,ist) = v1;
    SCM_v_2(:,ist) = v2;
    
    % Vertical velocity
    SCM_w(:,ist) = wmean;
    SCM_w_1(:,ist) = w1;
    SCM_w_2(:,ist) = w2;
    
    % Water: total, vapour, and liquid
    SCM_q(:,ist) = qmean;
    SCM_q_1(:,ist) = q1;
    SCM_q_2(:,ist) = q2;
    SCM_qv(:,ist) = qvmean;
    SCM_qv_1(:,ist) = vapour1;
    SCM_qv_2(:,ist) = vapour2;
    SCM_ql(:,ist) = qlmean;
    SCM_ql_1(:,ist) = liquid1;
    SCM_ql_2(:,ist) = liquid2;
    
    % (Liquid water) potential temperature
    SCM_th(:,ist) = thetamean;
    SCM_th_1(:,ist) = eos.theta1;
    SCM_th_2(:,ist) = eos.theta2;

    % Updraft buoyancy
    SCM_buoy(:,ist) = work.buoy;
    
    % Higher moments
    
    % TKE
    SCM_e_tot(:,ist) = tketot;
    SCM_e_res(:,ist) = (m1.*Rtke1 + m2.*Rtke2)./(m1 + m2);
    SCM_e1_res(:,ist) = Rtke1;
    SCM_e2_res(:,ist) = Rtke2;
    SCM_e_sg(:,ist) = (m1.*tke1 + m2.*tke2)./(m1 + m2);
    SCM_e1_sg(:,ist) = state_new.fluid(1).tke;
    SCM_e2_sg(:,ist) = state_new.fluid(2).tke;
    
    % u variance
    SCM_uu_tot(:,ist) = uutot;
    SCM_uu_res(:,ist) = (m1.*Ruu1 + m2.*Ruu2)./(m1 + m2);
    SCM_uu_1(:,ist) = Ruu1;
    SCM_uu_2(:,ist) = Ruu2;
    % Estimate SG u variance assuming isotropic turbulence
    SCM_uu_sg(:,ist) = (m1.*uu1 + m2.*uu2)./(m1 + m2);
    SCM_uu_sg1(:,ist) = uu1;
    SCM_uu_sg2(:,ist) = uu2;
    
    % v variance
    SCM_vv_tot(:,ist) = vvtot;
    SCM_vv_res(:,ist) = (m1.*Rvv1 + m2.*Rvv2)./(m1 + m2);
    SCM_vv_1(:,ist) = Rvv1;
    SCM_vv_2(:,ist) = Rvv2;
    % Estimate SG v variance assuming isotropic turbulence
    SCM_vv_sg(:,ist) = (m1.*vv1 + m2.*vv2)./(m1 + m2);
    SCM_vv_sg1(:,ist) = vv1;
    SCM_vv_sg2(:,ist) = vv2;
    
    % w variance
    SCM_ww_tot(:,ist) = wwtot;
    SCM_ww_res(:,ist) = (m1bar.*Rww1 + m2bar.*Rww2)./(m1bar + m2bar);
    SCM_ww_1(:,ist) = Rww1;
    SCM_ww_2(:,ist) = Rww2;
    % Estimate SG w variance assuming isotropic turbulence
    SCM_ww_sg(:,ist) = (m1.*ww1 + m2.*ww2)./(m1 + m2);
    SCM_ww_sg1(:,ist) = ww1;
    SCM_ww_sg2(:,ist) = ww2;
    
    % Water variance
    qq1 = state_new.fluid(1).varq;
    qq2 = state_new.fluid(2).varq;
    SCM_qq_tot(:,ist) = (m1.*(grid.aboves.*Rqq1(2:nzp) + grid.belows.*Rqq1(1:nz) + qq1) + m2.*(grid.aboves.*Rqq2(2:nzp) + grid.belows.*Rqq2(1:nz) + qq2))./(m1 + m2);;
    SCM_qq_res(:,ist) = (m1bar.*Rqq1 + m2bar.*Rqq2)./(m1bar + m2bar);
    SCM_qq_1(:,ist) = Rqq1;
    SCM_qq_2(:,ist) = Rqq2;
    SCM_qq_sg(:,ist) = (m1.*qq1 + m2.*qq2)./(m1 + m2);
    SCM_qq_sg1(:,ist) = qq1;
    SCM_qq_sg2(:,ist) = qq2;
    
    % Potential temperature variance
    thth1 = eos.Vartheta1;
    thth2 = eos.Vartheta2;
    SCM_thth_tot(:,ist) = (m1bar.*(Rthth1 + thth1) + m2bar.*(Rthth2 + thth2))./(m1bar + m2bar);;
    SCM_thth_res(:,ist) = (m1bar.*Rthth1 + m2bar.*Rthth2)./(m1bar + m2bar);
    SCM_thth_1(:,ist) = Rthth1;
    SCM_thth_2(:,ist) = Rthth2;
    % Approximate conversion from eta variance to theta variance
    % SCM_thth_sg1(:,ist) = state_new.fluid(1).vareta.*(eos.theta1p/Cpd).^2;
    % SCM_thth_sg2(:,ist) = state_new.fluid(2).vareta.*(eos.theta2p/Cpd).^2;
    SCM_thth_sg(:,ist) = (m1bar.*thth1 + m2bar.*thth2)./(m1bar + m2bar);
    SCM_thth_sg1(:,ist) = thth1;
    SCM_thth_sg2(:,ist) = thth2;
    
    % Moisture flux
    wq_res1 = (w1 - wmean).*(q1 - qmean);
    wq_res2 = (w2 - wmean).*(q2 - qmean);
    wq_sg1 = work.Dq1ed + work.Dq1bc;
    wq_sg2 = work.Dq2ed + work.Dq2bc;
    SCM_wq_res(:,ist) = (m1bar.*wq_res1 + m2bar.*wq_res2)./(m1bar + m2bar);
    SCM_wq_res1(:,ist) = wq_res1;
    SCM_wq_res2(:,ist) = wq_res2;
    SCM_wq_sg(:,ist) = (m1.*wq_sg1 + m2.*wq_sg2)./(m1 + m2);
    SCM_wq_sg1(:,ist) = wq_sg1;
    SCM_wq_sg2(:,ist) = wq_sg2;
    SCM_wq_tot(:,ist) = (m1bar.*(wq_res1 + weight_to_w(grid, wq_sg1)) + m2bar.*(wq_res2 + weight_to_w(grid, wq_sg2)))./(m1bar + m2bar);;
    
    % Potential temperature flux
    wth_res1 = (w1 - wmean).*(eos.theta1 - thetamean);
    wth_res2 = (w2 - wmean).*(eos.theta2 - thetamean);
    wth_sg1 = 0*thth1(2:nzp);
    wth_sg2 = 0*thth2(2:nzp);
    SCM_wth_res(:,ist) = (m1bar.*wth_res1 + m2bar.*wth_res2)./(m1bar + m2bar);
    SCM_wth_res1(:,ist) = wth_res1;
    SCM_wth_res2(:,ist) = wth_res2;
    SCM_wth_sg(:,ist) = (m1.*wth_sg1 + m2.*wth_sg2)./(m1 + m2);
    SCM_wth_sg1(:,ist) = wth_sg1;
    SCM_wth_sg2(:,ist) = wth_sg2;
    SCM_wth_tot(:,ist) = (m1bar.*(wth_res1 + weight_to_w(grid, wth_sg1)) + m2bar.*(wth_res2 + weight_to_w(grid, wth_sg2)))./(m1bar + m2bar);;
    
    % Entrainment
    SCM_M21_instab(:,ist) = relabel.M21_instab;
    SCM_M21_sort(:,ist) = relabel.M21_sort;
    SCM_M21_dwdz(:,ist) = relabel.M21_dwdz;
    SCM_M21_mix(:,ist) = relabel.M21_mix;
    SCM_M21(:,ist) = relabel.M21;
    
    % Detrainment
    SCM_M12_instab(:,ist) = relabel.M12_instab;
    SCM_M12_sort(:,ist) = relabel.M12_sort;
    SCM_M12_dwdz(:,ist) = relabel.M12_dwdz;
    SCM_M12_mix(:,ist) = relabel.M12_mix;
    SCM_M12(:,ist) = relabel.M12;
    
end


% At the end of the run write to file
if istep == time.nstop

    SCM_zw = grid.zw;
    SCM_zp = grid.zp;
    
    SCM_time_ser = ts.time_high_res;
    SCM_zstar = ts.zstar;
    SCM_zcbase = ts.zcbaseSG;
    SCM_zctop = ts.zctopSG;
    SCM_cldcov = ts.totcldcov;
    SCM_cloud_fraction  = ts.cloud_fraction;
    SCM_cloud_fraction1 = ts.cloud_fraction1;
    SCM_cloud_fraction2 = ts.cloud_fraction2;
    SCM_cloud_fraction1_sigma1 = ts.cloud_fraction1_sigma1;
    SCM_cloud_fraction2_sigma2 = ts.cloud_fraction2_sigma2;

    % Write to file
    filename = fullfile(settings.folders.data_scm, '2FSCM_results.mat');
    % Save all variables whose name begins SCM...
    save(filename,'-regexp','^SCM');
    
    % Remove the "SCM_" from the variable names
    structure_old = load(filename);
    fieldnames_old = fieldnames(structure_old);
    fieldnames_new = strrep(fieldnames_old, 'SCM_', ''); 
    for z = 1:length(fieldnames_old)
      structure_new.(fieldnames_new{z}) = structure_old.(fieldnames_old{z});
    end
    save(filename,'-struct','structure_new');
end


