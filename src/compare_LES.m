% Save data for direct comparison with LES results

% List of times at which to compare with LES
SCM_times = [14000, 21800, 32600, 42800];
nst = numel(SCM_times); 

% Time step size
dt = time.dt;


% Is it time to save data?
lsave = 0;
ist = 0;
while ist < nst & ~lsave
    ist = ist + 1;
    st = SCM_times(ist);
    if (time.t - dt < st & time.t >= st)
        lsave = 1;
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
    Rww1 = (wmean - w1).^2;
    Rww2 = (wmean - w2).^2;
    % Resolved TKE per unit mass
    Rtke1 = 0.5*((umean - u1).^2 + (vmean - v1).^2 + grid.aboves.*Rww1(2:nzp) + grid.belows.*Rww1(1:nz));
    Rtke2 = 0.5*((umean - u2).^2 + (vmean - v2).^2 + grid.aboves.*Rww2(2:nzp) + grid.belows.*Rww2(1:nz));
    % Total TKE per unit mass
    tketot = (m1.*(Rtke1 + tke1) + m2.*(Rtke2 + tke2))./(m1 + m2);
    
    % Resolved q variance
    qmean = (m1bar.*q1 + m2bar.*q2)./(m1bar + m2bar);
    Rqq1 = (q1 - qmean).^2;
    Rqq2 = (q2 - qmean).^2;
    
    % Resolved theta(_l) variance.
    thetamean = (m1bar.*eos.theta1 + m2bar.*eos.theta2)./(m1bar + m2bar);
    Rthth1 = (eos.theta1 - thetamean).^2;
    Rthth2 = (eos.theta2 - thetamean).^2;
    
    % Pack data ready for output
    
    % Mean profiles
    
    % Mass fractions at p levels and w levels
    SCM_sigma1(:,ist) = sigma1;
    SCM_sigma2(:,ist) = sigma2;
    SCM_sigma1w(:,ist) = sigma1w;
    SCM_sigma2w(:,ist) = sigma2w;
    
    % Vertical velocity
    SCM_w_1(:,ist) = w1;
    SCM_w_2(:,ist) = w2;
    
    % Water: total, vapour, and liquid
    SCM_q_1(:,ist) = q1;
    SCM_q_2(:,ist) = q2;
    SCM_qv_1(:,ist) = vapour1;
    SCM_qv_2(:,ist) = vapour2;
    SCM_ql_1(:,ist) = liquid1;
    SCM_ql_2(:,ist) = liquid2;
    
    % (Liquid water) potential temperature
    SCM_th_1(:,ist) = eos.theta1;
    SCM_th_2(:,ist) = eos.theta2;

    % Updraft buoyancy
    SCM_buoy(:,ist) = work.buoy;
    
    % Higher moments
    
    % TKE
    SCM_e1_res(:,ist) = Rtke1;
    SCM_e2_res(:,ist) = Rtke2;
    SCM_e1_sg(:,ist) = state_new.fluid(1).tke;
    SCM_e2_sg(:,ist) = state_new.fluid(2).tke;
    
    % w variance
    SCM_ww1(:,ist) = Rww1;
    SCM_ww2(:,ist) = Rww2;
    % Estimate SG w variance assuming isotropic turbulence
    SCM_ww_sg1(:,ist) = (2/3)*state_new.fluid(1).tke;
    SCM_ww_sg2(:,ist) = (2/3)*state_new.fluid(2).tke;
    
    % Water variance
    SCM_qq1(:,ist) = Rqq1;
    SCM_qq2(:,ist) = Rqq2;
    SCM_qq_sg1(:,ist) = state_new.fluid(1).varq;
    SCM_qq_sg2(:,ist) = state_new.fluid(2).varq;
    
    % Potential temperature variance
    SCM_thth1(:,ist) = Rthth1;
    SCM_thth2(:,ist) = Rthth2;
    % Approximate conversion from eta variance to theta variance
    SCM_thth_sg1(:,ist) = state_new.fluid(1).vareta.*(eos.theta1/Cpd).^2;
    SCM_thth_sg2(:,ist) = state_new.fluid(2).vareta.*(eos.theta2/Cpd).^2;
    
end


% At the end of the run write to file
if istep == time.nstop

    SCM_zw = grid.zw;
    SCM_zp = grid.zp;
    
    SCM_time_ser = ts.time;
    if exist('zcldtop')
        % This will be true if plot_time_series has been called
        SCM_zstar = ts.zstar;
        SCM_zcbase = ts.zcbaseSG;
        SCM_zctop = ts.zctopSG;
        SCM_cldcov = ts.totcldcov;
    end

    % Write to file
    filename = fullfile(settings.folders.data_scm, 'SCM_results.mat');
    % Save all variables whose name begins SCM...
    save(filename,'-regexp','^SCM');

end


