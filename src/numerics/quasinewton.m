% Take a number of quasi-Newton iteration


% Levels for displaying diagnostics
%krange = nz-4:nz;
krange = 56:58; %53:58;
% and for plotting
prange = 20:50;

% Unpack some fields that are needed several times
nz = grid.nz;
nzp = nz + 1;
dzp = grid.dzp;
dzw = grid.dzw;
abovep = grid.abovep;
belowp = grid.belowp;
abovew = grid.abovew;
beloww = grid.beloww;
abover = grid.abover;
belowr = grid.belowr;
aboves = grid.aboves;
belows = grid.belows;
% Useful quantities for picking the middle coefficient for
% relabelling terms
arbp = abover(1:nz ).*belowp;
brap = belowr(2:nzp).*abovep;

% Fixer to ensure tke does not drop below specified minimum.
% It is initialized to zero but may be modified as iterations proceed
fixmtke1 = zeros(1,nz);
fixmtke2 = zeros(1,nz);

% ------

for qn_iter = 1:qn_iter_max

    % disp(' ')
    disp(['Iteration ' num2str(qn_iter)])
 
    % Check for nans and infs in state_new
    %check_this = 'state';
    %check_for_nans
    
    % ------
    
    % Investigate sensitivity of residuals to perturbing unknowns
%     if istep == 3
%         sensitivity
%     end

    % ------
    
    if conv_diag
        
        % For convergence testing, rescale error
        rescale_error
  
        % Pack errors into one array
        xx(1:9:9*nz+1) = state_err.fluid(1).w;
        xx(2:9:9*nz+2) = state_err.fluid(2).w;
        xx(3:9:9*nz+3) = state_err.fluid(1).eta;
        xx(4:9:9*nz+4) = state_err.fluid(2).eta;
        xx(5:9:9*nz+5) = state_err.fluid(1).q;
        xx(6:9:9*nz+6) = state_err.fluid(2).q;
        xx(7:9:9*nz-2) = state_err.fluid(1).m;
        xx(8:9:9*nz-1) = state_err.fluid(2).m;
        xx(9:9:9*nz)   = state_err.p;
        
        % Consider doing the same for tke and variance errors
        % Need to compute and scale variance errors in rescale_error
        
    end
    
    
    % Tendencies computed from new state
    dt = time.dt;
    t_new = time.t + dt;
    old_diff.flag = 0;
    [tend,relabel,eos,force,scales,surface_flux,budgets,work] = ...
              tendencies(grid,state_new,settings,t_new,dt,switches,old_diff);
    
    % For debugging
    % if istep > 17
    %  dump_tend
    % end
    
    % Check for nans and infs in tendencies
    %check_this = 'tend';
    %check_for_nans
    
    % Plot variance budgets; they should balance at convergence
    % plot_var_budgets
    
    % Overwrite selected tendencies for checking convergence
    % overwrite_tendencies
    

    % Unpack some fields for clarity of code
    m1 = state_new.fluid(1).m;
    m2 = state_new.fluid(2).m;
    w1 = state_new.fluid(1).w;
    w2 = state_new.fluid(2).w;
    eta1 = state_new.fluid(1).eta;
    eta2 = state_new.fluid(2).eta;
    q1 = state_new.fluid(1).q;
    q2 = state_new.fluid(2).q;
    u1 = state_new.fluid(1).u;
    u2 = state_new.fluid(2).u;
    v1 = state_new.fluid(1).v;
    v2 = state_new.fluid(2).v;
    tke1 = state_new.fluid(1).tke;
    tke2 = state_new.fluid(2).tke;
    m1bar = work.m1bar;
    m2bar = work.m2bar;
    dpdz = work.dpdz;
    M12 = relabel.M12;
    M21 = relabel.M21;
    M12bar = relabel.M12bar;
    M21bar = relabel.M21bar;
    what12_w1     = relabel.what12   - work.w1ubar;
    what21_w1     = relabel.what21   - work.w1ubar;
    what12_w2     = relabel.what12   - work.w2ubar;
    what21_w2     = relabel.what21   - work.w2ubar;
    etahat12_eta1 = relabel.etahat12 - work.eta1ubar;
    etahat21_eta1 = relabel.etahat21 - work.eta1ubar;
    etahat12_eta2 = relabel.etahat12 - work.eta2ubar;
    etahat21_eta2 = relabel.etahat21 - work.eta2ubar;
    qhat12_q1     = relabel.qhat12   - work.q1ubar;
    qhat21_q1     = relabel.qhat21   - work.q1ubar;
    qhat12_q2     = relabel.qhat12   - work.q2ubar;
    qhat21_q2     = relabel.qhat21   - work.q2ubar;

    
    % Left hand sides of all equations
    adt = time.alpha*time.dt;
    if switches.c
        lhs1m   = m2                            + m1                            - adt*tend.fluid(1).m.tot;
        lhs1eta = state_new.fluid(2).eta.*m2bar + state_new.fluid(1).eta.*m1bar - adt*tend.fluid(1).meta.tot;
        lhs1q   = state_new.fluid(2).q  .*m2bar + state_new.fluid(1).q  .*m1bar - adt*tend.fluid(1).mq.tot;
        lhs1w   = state_new.fluid(2).w  .*m2bar + state_new.fluid(1).w  .*m1bar - adt*tend.fluid(1).mw.tot;
        lhs1u   = state_new.fluid(2).u  .*m2    + state_new.fluid(1).u  .*m1    - adt*tend.fluid(1).mu.tot;
        lhs1v   = state_new.fluid(2).v  .*m2    + state_new.fluid(1).v  .*m1    - adt*tend.fluid(1).mv.tot;
        lhs1tke = state_new.fluid(2).tke.*m2    + state_new.fluid(1).tke.*m1    - adt*tend.fluid(1).mtke.tot;
        lhs2m   = - adt*tend.fluid(2).m.tot;
        lhs2eta = - adt*tend.fluid(2).meta.tot;
        lhs2q   = - adt*tend.fluid(2).mq.tot;
        lhs2w   = - adt*tend.fluid(2).mw.tot;
        lhs2u   = - adt*tend.fluid(2).mu.tot;
        lhs2v   = - adt*tend.fluid(2).mv.tot;
        lhs2tke = - adt*tend.fluid(2).mtke.tot;
        disp('** include tke fixer **')
        pause
    else
        lhs1m   = m1                            - adt*tend.fluid(1).m.tot;
        lhs1eta = state_new.fluid(1).eta.*m1bar - adt*tend.fluid(1).meta.tot;
        lhs1q   = state_new.fluid(1).q  .*m1bar - adt*tend.fluid(1).mq.tot;
        lhs1w   = state_new.fluid(1).w  .*m1bar - adt*tend.fluid(1).mw.tot;
        lhs1u   = state_new.fluid(1).u  .*m1    - adt*tend.fluid(1).mu.tot;
        lhs1v   = state_new.fluid(1).v  .*m1    - adt*tend.fluid(1).mv.tot;
        lhs1tke = state_new.fluid(1).tke.*m1    - adt*tend.fluid(1).mtke.tot - fixmtke1;
        lhs2m   = m2                            - adt*tend.fluid(2).m.tot;
        lhs2eta = state_new.fluid(2).eta.*m2bar - adt*tend.fluid(2).meta.tot;
        lhs2q   = state_new.fluid(2).q  .*m2bar - adt*tend.fluid(2).mq.tot;
        lhs2w   = state_new.fluid(2).w  .*m2bar - adt*tend.fluid(2).mw.tot;
        lhs2u   = state_new.fluid(2).u  .*m2    - adt*tend.fluid(2).mu.tot;
        lhs2v   = state_new.fluid(2).v  .*m2    - adt*tend.fluid(2).mv.tot;
        lhs2tke = state_new.fluid(2).tke.*m2    - adt*tend.fluid(2).mtke.tot - fixmtke2;
    end
    
    
    % Residuals in all equations (note opposite sign convention
    % to Gibbs paper)
    res1m   =  rhs1m   - lhs1m;
    res1eta =  rhs1eta - lhs1eta;
    res1q   =  rhs1q   - lhs1q;
    res1w   =  rhs1w   - lhs1w;
    res1u   =  rhs1u   - lhs1u;
    res1v   =  rhs1v   - lhs1v;
    res1tke =  rhs1tke - lhs1tke;
    res2m   =  rhs2m   - lhs2m;
    res2eta =  rhs2eta - lhs2eta;
    res2q   =  rhs2q   - lhs2q;
    res2w   =  rhs2w   - lhs2w;
    res2u   =  rhs2u   - lhs2u;
    res2v   =  rhs2v   - lhs2v;
    res2tke =  rhs2tke - lhs2tke;

    
    if conv_diag
        
        % Fix residuals to mimic converged solution
        res1m   = res1m   + resfix1m;
        res1eta = res1eta + resfix1eta;
        res1q   = res1q   + resfix1q;
        res1w   = res1w   + resfix1w;
        res1u   = res1u   + resfix1u;
        res1v   = res1v   + resfix1v;
        res1tke = res1tke + resfix1tke;
        res2m   = res2m   + resfix2m;
        res2eta = res2eta + resfix2eta;
        res2q   = res2q   + resfix2q;
        res2w   = res2w   + resfix2w;
        res2u   = res2u   + resfix2u;
        res2v   = res2v   + resfix2v;
        res2tke = res2tke + resfix2tke;
        eos.res_eta1 = eos.res_eta1 + resfixeoseta1;
        eos.res_eta2 = eos.res_eta2 + resfixeoseta2;
        eos.res_etap1 = eos.res_etap1 + resfixeosetap1;
        eos.res_etap2 = eos.res_etap2 + resfixeosetap2;
        eos.res_rho1 = eos.res_rho1 + resfixeosrho1;
        eos.res_rho2 = eos.res_rho2 + resfixeosrho2;
        eos.res_sigma = eos.res_sigma + resfixeossigma;
        
        % Disect residuals into contributions from different terms
        % disect

    end
    
    
    % Interpolate mass residuals to w levels
    % using `reversed' weighting for conservation
    res1mbar = weight_to_w(grid,res1m);
    res2mbar = weight_to_w(grid,res2m);
    
    
    % Residuals in w equation after substituting for buoyancy
    res1w = res1w + adt*m1bar.*dpdz.*eos.drdeta1.*eos.res_eta1;
    res2w = res2w + adt*m2bar.*dpdz.*eos.drdeta2.*eos.res_eta2;
    
    % Implied residuals in quasi-advective form eta, q, w, u and v equations
    res1eta = res1eta - work.eta1ubar.*res1mbar;
    res2eta = res2eta - work.eta2ubar.*res2mbar;
    res1q   = res1q   - work.q1ubar  .*res1mbar;
    res2q   = res2q   - work.q2ubar  .*res2mbar;
    res1w   = res1w   - work.w1ubar  .*res1mbar;
    res2w   = res2w   - work.w2ubar  .*res2mbar;
    res1u   = res1u   - work.u1ubar  .*res1m;
    res2u   = res2u   - work.u2ubar  .*res2m;
    res1v   = res1v   - work.v1ubar  .*res1m;
    res2v   = res2v   - work.v2ubar  .*res2m;
    
    % Residual in sigma equation after substituting for temperature
    res_s = eos.res_sigma + m1.*(eos.res_rho1 + eos.drdetap1.*eos.res_etap1) ...
                          + m2.*(eos.res_rho2 + eos.drdetap2.*eos.res_etap2);
                      
    
    % *** for testing ***
    %disp('*** modified residuals ***')
    %res1m = 0*res1m;
    %res2m = 0*res2m;
    %res1eta = 0*res1eta;
    %res2eta = 0*res2eta;
    %res1q = 0*res1q;
    %res2q = 0*res2q;
    %res1w = 0*res1w;
    %res2w = 0*res2w;
                      
%    disp(' ')
%    disp('Residuals')
%     disp(['res_w1    ' num2str(res1w(krange))])
%     disp(['res_w2    ' num2str(res2w(krange))])
%     disp(['res_m1    ' num2str(res1m(krange))])
%     disp(['res_m2    ' num2str(res2m(krange))])
%     disp(['res_eta1  ' num2str(res1eta(krange))])
%     disp(['res_eta2  ' num2str(res2eta(krange))])
%     disp(['res_q1    ' num2str(res1q(krange))])
%     disp(['res_q2    ' num2str(res2q(krange))])
%     disp(['res_u1    ' num2str(res1u(krange))])
%     disp(['res_u2    ' num2str(res2u(krange))])
%     disp(['res_v1    ' num2str(res1v(krange))])
%     disp(['res_v2    ' num2str(res2v(krange))])
%     disp(['res_sigma ' num2str(res_s(krange))])
%     disp(['res_etap1 ' num2str(eos.res_etap1(krange))])
%     disp(['res_eta1  ' num2str(eos.res_eta1(krange))])
%     disp(['res_etap2 ' num2str(eos.res_etap2(krange))])
%     disp(['res_eta2  ' num2str(eos.res_eta2(krange))])
%    disp(['res_tke1  ' num2str(res1tke(krange))])
%    disp(['res_tke2  ' num2str(res2tke(krange))])
%     disp(['res_Veta1  ' num2str(tend.fluid(1).mvareta.tot(krange))])
%     disp(['res_Veta2  ' num2str(tend.fluid(2).mvareta.tot(krange))])
%     disp(['res_Vq1    ' num2str(tend.fluid(1).mvarq.tot(krange))])
%     disp(['res_Vq2    ' num2str(tend.fluid(2).mvarq.tot(krange))])

%     disp('Max eta var tend')
%     [tx,kx] = max(tend.fluid(1).mvareta.tot)
%     [tx,kx] = max(tend.fluid(2).mvareta.tot)     
%     disp('Max q var tend')
%     [tx,kx] = max(tend.fluid(1).mvarq.tot)
%     [tx,kx] = max(tend.fluid(2).mvarq.tot)
%    disp('Max tke res')
%    [tx,kx] = max(res1tke)
%    [tx,kx] = max(res2tke)
%    disp('Max m res')
%    [tx,kx] = max(res1m)
%    [tx,kx] = max(res2m)

% Plots of max residuals vs iteration for checking convergence
save_res_convergence
% if (qn_iter == qn_iter_max)
%     plot_res_convergence
%     pause
% end

    if conv_diag
        fs = 15;
        figure(10)
        subplot(2,3,1)
        plot(res_s(prange),grid.zp(prange))
        title('s res')
        set(gca,'fontsize',fs)
        subplot(2,3,2)
        plot(res1w(prange),grid.zw(prange),'b',res2w(prange),grid.zw(prange),'r')
        title('w res')
        set(gca,'fontsize',fs)
        subplot(2,3,3)
        plot(res1m(prange),grid.zp(prange),'b',res2m(prange),grid.zp(prange),'r')
        title('m res')
        set(gca,'fontsize',fs)
        subplot(2,3,4)
        plot(res1eta(prange),grid.zw(prange),'b',res2eta(prange),grid.zw(prange),'r')
        title('eta res')
        set(gca,'fontsize',fs)
        subplot(2,3,5)
        plot(res1q(prange),grid.zw(prange),'b',res2q(prange),grid.zw(prange),'r')
        title('q res')
        set(gca,'fontsize',fs)
        %figure(1)
    end
    
    
    % --------
    
    % Now build the matrix   
    build_linear_system
    
    % --------
    
    if conv_diag
        
        % For testing, compute predicted residuals using linearization
        % to compare with actual residuals ...
        rr = - Ndiagmult(cc,xx);
        rr_w1   = rr(1:9:9*nz+1);
        rr_w2   = rr(2:9:9*nz+2);
        rr_eta1 = rr(3:9:9*nz+3);
        rr_eta2 = rr(4:9:9*nz+4);
        rr_q1   = rr(5:9:9*nz+5);
        rr_q2   = rr(6:9:9*nz+6);
        rr_m1   = rr(7:9:9*nz-2);
        rr_m2   = rr(8:9:9*nz-1);
        rr_s    = rr(9:9:9*nz);

        disp(' ')
        disp('Predicted Residuals')
        disp(['res_w1    ' num2str(rr_w1(krange))])
        disp(['res_w2    ' num2str(rr_w2(krange))])
        disp(['res_m1    ' num2str(rr_m1(krange))])
        disp(['res_m2    ' num2str(rr_m2(krange))])
        disp(['res_eta1  ' num2str(rr_eta1(krange))])
        disp(['res_eta2  ' num2str(rr_eta2(krange))])
        disp(['res_q1    ' num2str(rr_q1(krange))])
        disp(['res_q2    ' num2str(rr_q2(krange))])
        disp(['res_sigma ' num2str(rr_s(krange))])
    
        figure(11)
        subplot(2,3,1)
        plot(rr_s(prange),grid.zp(prange))
        title('pred. s res')
        set(gca,'fontsize',fs)
        subplot(2,3,2)
        plot(rr_w1(prange),grid.zw(prange),'b',rr_w2(prange),grid.zw(prange),'r')
        title('pred. w res')
        set(gca,'fontsize',fs)
        subplot(2,3,3)
        plot(rr_m1(prange),grid.zp(prange),'b',rr_m2(prange),grid.zp(prange),'r')
        title('pred. m res')
        set(gca,'fontsize',fs)
        subplot(2,3,4)
        plot(rr_eta1(prange),grid.zw(prange),'b',rr_eta2(prange),grid.zw(prange),'r')
        title('pred. eta res')
        set(gca,'fontsize',fs)
        subplot(2,3,5)
        plot(rr_q1(prange),grid.zw(prange),'b',rr_q2(prange),grid.zw(prange),'r')
        title('pred. q res')
        set(gca,'fontsize',fs)
        %figure(1)
        
        % Disect predicted residuals
        % disect_prediction
    
        % For testing: replace true residual by predicted residual
%         disp('** residuals replaced **')
%         res1w = rr_w1;
%         res2w = rr_w2;
%         res1eta = rr_eta1;
%         res2eta = rr_eta2;
%         res1q = rr_q1;
%         res2q = rr_q2;
%         res1m = rr_m1;
%         res2m = rr_m2;
%         res_s = rr_s;
        
    end
    
    % --------
    
    % Build linear system for w-eta-q-m-p system
    
    % Now build the right hand side
    ix = 1:9:9*nz+1;
    rhs(ix) = res1w;
    ix = ix + 1;
    rhs(ix) = res2w;
    ix = ix + 1;
    rhs(ix) = res1eta;
    ix = ix + 1;
    rhs(ix) = res2eta;
    ix = ix + 1;
    rhs(ix) = res1q;
    ix = ix + 1;
    rhs(ix) = res2q;
    ix = 7:9:9*(nz-1)+7;
    rhs(ix) = res1m;
    ix = ix + 1;
    rhs(ix) = res2m;
    ix = ix + 1;
    rhs(ix) = res_s;
    
    % ---------
    
    % Look at which terms in cc inverse contribute to
    % certain increments
    %invert_cc(cc,rhs)
    %pause

    % ---------
    
    % Now solve the linear system
    
    % FindEvec(cc, grid)
    
%     % check for nans and infs
%     nan_in_cc = sum(sum(isnan(cc)));
%     nan_in_rhs = sum(isnan(rhs));
%     inf_in_cc = sum(sum(isinf(cc)));
%     inf_in_rhs = sum(isinf(rhs));
%     if sum([nan_in_cc,nan_in_rhs,inf_in_cc,inf_in_rhs])
%         disp(['nan_in_cc  = ',num2str(nan_in_cc)])
%         disp(['nan_in_rhs = ',num2str(nan_in_rhs)])
%         disp(['inf_in_cc  = ',num2str(inf_in_cc)])
%         disp(['inf_in_rhs = ',num2str(inf_in_rhs)])
%         pause
%     end
    
    xx = Ndiagsolveb(cc,rhs);
    %xx = Ndiagsolvex(cc,rhs);
    
    % ---------
    
    % and unpack the increments
    
    % disp('*** half eta increments ***')
    inc_w1   = xx(1:9:9*nz+1);
    inc_w2   = xx(2:9:9*nz+2);
    inc_eta1 = xx(3:9:9*nz+3);
    inc_eta2 = xx(4:9:9*nz+4);
    inc_q1   = xx(5:9:9*nz+5);
    inc_q2   = xx(6:9:9*nz+6);
    inc_m1   = xx(7:9:9*nz-2);
    inc_m2   = xx(8:9:9*nz-1);
    inc_p    = xx(9:9:9*nz);
    
    % ---------
    
%     % Bound mass increments to prevent negative mass, preserve conservation
     corrm1 = max(0,-0.9*m1 - inc_m1);
     corrm2 = max(0,-0.9*m2 - inc_m2);
     inc_m1 = inc_m1 + corrm1 - corrm2;
     inc_m2 = inc_m2 + corrm2 - corrm1;
     
     if sum(corrm1 + corrm2)
         disp('** Mass correction made **')
         % pause
     end
    
    % ---------
    
    % Backsubstitute for T increments
    
    % T increments at p levels
    inc_eta_bar = aboves.*inc_eta1(2:nzp) + belows.*inc_eta1(1:nz);
    inc_q_bar   = aboves.*inc_q1(2:nzp)   + belows.*inc_q1(1:nz);
    for k = 1:nz
        qbar = aboves(k)*state_new.fluid(1).q(k+1) + belows(k)*state_new.fluid(1).q(k);
        [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = ...
                    gibbs(state_new.p(k),state_new.fluid(1).T(k),qbar,constants.therm);
        inc_T1(k) = -(eos.res_etap1(k) + inc_eta_bar(k) + gpt*inc_p(k) + gtw*inc_q_bar(k))/gtt;
    end
    inc_eta_bar = aboves.*inc_eta2(2:nzp) + belows.*inc_eta2(1:nz);
    inc_q_bar   = aboves.*inc_q2(2:nzp)   + belows.*inc_q2(1:nz);
    for k = 1:nz
        qbar = aboves(k)*state_new.fluid(2).q(k+1) + belows(k)*state_new.fluid(2).q(k);
        [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = ...
                    gibbs(state_new.p(k),state_new.fluid(2).T(k),qbar,constants.therm);
        inc_T2(k) = -(eos.res_etap2(k) + inc_eta_bar(k) + gpt*inc_p(k) + gtw*inc_q_bar(k))/gtt;
    end
    
    % T increments at w levels
    for k = 1:nzp
        
        if k == 1
            pbar      = grid.extrapb1*state_new.p(1) + grid.extrapb2*state_new.p(2);
            inc_p_bar = grid.extrapb1*inc_p(1)   + grid.extrapb2*inc_p(2);
        elseif k == nzp
            pbar      = grid.extraptnz*state_new.p(nz) + grid.extraptnzm*state_new.p(nz-1);
            inc_p_bar = grid.extraptnz*inc_p(nz)   + grid.extraptnzm*inc_p(nz-1);
        else
            pbar      = grid.abovew(k)*state_new.p(k) ...
                      + grid.beloww(k)*state_new.p(k-1);
            inc_p_bar = grid.abovew(k)*inc_p(k) ...
                      + grid.beloww(k)*inc_p(k-1);
        end
        [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = ...
                    gibbs(pbar,state_new.fluid(1).Tw(k),state_new.fluid(1).q(k),constants.therm);
        inc_Tw1(k) = -(eos.res_eta1(k) + inc_eta1(k) + gpt*inc_p_bar + gtw*inc_q1(k))/gtt;
        [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = ...
                    gibbs(pbar,state_new.fluid(2).Tw(k),state_new.fluid(2).q(k),constants.therm);
        inc_Tw2(k) = -(eos.res_eta2(k) + inc_eta2(k) + gpt*inc_p_bar + gtw*inc_q2(k))/gtt;
 
    end
    
    
%     disp([' '])
%     disp('Increments')
%     disp(['inc_p     ' num2str(inc_p(krange))])
%     disp(['inc_w1    ' num2str(inc_w1(krange))])
%     disp(['inc_w2    ' num2str(inc_w2(krange))])
%     disp(['inc_w2    ' num2str(inc_w2(52))])
%     disp(['inc_m1    ' num2str(inc_m1(krange))])
%     disp(['inc_m2    ' num2str(inc_m2(krange))])
%     disp(['inc_eta1  ' num2str(inc_eta1(krange))])
%     disp(['inc_eta2  ' num2str(inc_eta2(krange))])
%     disp(['inc_eta2  ' num2str(inc_eta2(52))])
%     disp(['inc_q1    ' num2str(inc_q1(krange))])
%     disp(['inc_q2    ' num2str(inc_q2(krange))])
%     disp(['inc_T1    ' num2str(inc_T1(krange))])
%     disp(['inc_T2    ' num2str(inc_T2(krange))])
%     disp(['inc_Tw1   ' num2str(inc_Tw1(krange))])
%     disp(['inc_Tw2   ' num2str(inc_Tw2(krange))])
    
    if conv_diag
        disp([' '])
        disp('Increments')
        disp(['inc_p     ' num2str(inc_p(krange))])
        disp(['inc_w1    ' num2str(inc_w1(krange))])
        disp(['inc_w2    ' num2str(inc_w2(krange))])
        disp(['inc_m1    ' num2str(inc_m1(krange))])
        disp(['inc_m2    ' num2str(inc_m2(krange))])
        disp(['inc_eta1  ' num2str(inc_eta1(krange))])
        disp(['inc_eta2  ' num2str(inc_eta2(krange))])
        disp(['inc_q1    ' num2str(inc_q1(krange))])
        disp(['inc_q2    ' num2str(inc_q2(krange))])
        disp(['inc_T1    ' num2str(inc_T1(krange))])
        disp(['inc_T2    ' num2str(inc_T2(krange))])
        disp(['inc_Tw1   ' num2str(inc_Tw1(krange))])
        disp(['inc_Tw2   ' num2str(inc_Tw2(krange))])
        
        figure(9)
        fs = 15;
        subplot(2,3,1)
        plot(inc_p(prange),grid.zp(prange))
        title('p inc')
        set(gca,'fontsize',fs)
        subplot(2,3,2)
        plot(inc_w1(prange),grid.zw(prange),'b',inc_w2(prange),grid.zw(prange),'r')
        title('w inc')
        set(gca,'fontsize',fs)
        subplot(2,3,3)
        plot(inc_m1(prange),grid.zp(prange),'b',inc_m2(prange),grid.zp(prange),'r')
        title('m inc')
        set(gca,'fontsize',fs)
        subplot(2,3,4)
        plot(inc_eta1(prange),grid.zw(prange),'b',inc_eta2(prange),grid.zw(prange),'r')
        title('eta inc')
        set(gca,'fontsize',fs)
        subplot(2,3,5)
        plot(inc_q1(prange),grid.zw(prange),'b',inc_q2(prange),grid.zw(prange),'r')
        title('q inc')
        set(gca,'fontsize',fs)
        subplot(2,3,6)
        plot(inc_T1(prange),grid.zp(prange),'b',inc_T2(prange),grid.zp(prange),'r',...
             inc_Tw1(prange),grid.zw(prange),'b--',inc_Tw2(prange),grid.zw(prange),'r--')
        title('T inc')
        set(gca,'fontsize',fs)
        pause
        figure(1)
        
    end

    
    % Increment all variables
    state_new.p            =     state_new.p            + inc_p;
    state_new.fluid(1).m   =     state_new.fluid(1).m   + inc_m1;
    state_new.fluid(1).w   =     state_new.fluid(1).w   + inc_w1;
    state_new.fluid(1).eta =     state_new.fluid(1).eta + inc_eta1;
    state_new.fluid(1).q   = max(state_new.fluid(1).q   + inc_q1, 0);
    state_new.fluid(1).T   =     state_new.fluid(1).T   + inc_T1;
    state_new.fluid(1).Tw  =     state_new.fluid(1).Tw  + inc_Tw1;
    state_new.fluid(2).m   =     state_new.fluid(2).m   + inc_m2;
    state_new.fluid(2).w   =     state_new.fluid(2).w   + inc_w2;
    state_new.fluid(2).eta =     state_new.fluid(2).eta + inc_eta2;
    state_new.fluid(2).q   = max(state_new.fluid(2).q   + inc_q2, 0);
    state_new.fluid(2).T   =     state_new.fluid(2).T   + inc_T2;
    state_new.fluid(2).Tw  =     state_new.fluid(2).Tw  + inc_Tw2;
    
    if min(inc_Tw2) < -0.1
        disp(['Min inc_Tw2 = ',num2str(min(inc_Tw2))])
        disp(['massterm2 = ',num2str(massterm2)])
        %pause
    end
    
    % --------
    
    % Build linear system for u1 system
    ix = 1:nz;
    rhsu = res1u;
    
    dd = zeros(3,nz);
    % Tendency term
    dd(2,ix) = dd(2,ix) + m1;
    % Transport terms    
    dd(1,ix) = dd(1,ix) - adt*work.dFu1dub(1:nz)./dzp;
    dd(2,ix) = dd(2,ix) - adt*(work.dFu1dua(1:nz) - work.dFu1dub(2:nzp))./dzp;
    dd(3,ix) = dd(3,ix) + adt*work.dFu1dua(2:nzp)./dzp;
    % Diffusion terms
    dd(1,ix) = dd(1,ix) - adt*work.dDu1dub(1:nz)./dzp;
    dd(2,ix) = dd(2,ix) - adt*(work.dDu1dua(1:nz) - work.dDu1dub(2:nzp))./dzp;
    dd(3,ix) = dd(3,ix) + adt*work.dDu1dua(2:nzp)./dzp;
    % Relabelling terms - only keep dependence on u1
    dd(2,ix) = dd(2,ix) - adt*M12.*relabel.duhat12du1 ...
                        + adt*M21.*relabel.duhat21du1;
    
    % Now solve the linear system
    % disp('*** zero increments ***')
    inc_u1 = Ndiagsolveb(dd,rhsu);
    %inc_u1 = Ndiagsolvex(dd,rhsu);
    
    if conv_diag
        % For testing, compute predicted residuals using linearization
        % to compare with actual residuals ...
        yy = state_err.fluid(1).u;
        rr = - Ndiagmult(dd,yy);
        disp(' ')
        disp('Predicted Residual')
        disp(['res_u1    ' num2str(rr(krange))])
        zz = - adt*M12.*relabel.duhat12du2 + adt*M21.*relabel.duhat21du2;
        if switches.c
            zz = zz + m2;
        end
        rrz = - zz.*state_err.fluid(2).u;
        rr = rr + rrz;
        disp(['with u2pert ',num2str(rr(krange))])
        %pause
    end
 
    % And increment u1
    state_new.fluid(1).u = state_new.fluid(1).u + inc_u1;
    if conv_diag
        disp(['inc_u1 = ',num2str(inc_u1(krange))])
    end
    
    % --------
    
    % Build linear system for u2 system
    ix = 1:nz;
    rhsu = res2u;
    
    dd = zeros(3,nz);
    % Tendency term
    if ~switches.c
        dd(2,ix) = dd(2,ix) + m2;
    else
        dd(2,ix) = dd(2,ix) + 0.01*m2;
    end
    % Transport terms    
    dd(1,ix) = dd(1,ix) - adt*work.dFu2dub(1:nz)./dzp;
    dd(2,ix) = dd(2,ix) - adt*(work.dFu2dua(1:nz) - work.dFu2dub(2:nzp))./dzp;
    dd(3,ix) = dd(3,ix) + adt*work.dFu2dua(2:nzp)./dzp;
    % Diffusion terms
    dd(1,ix) = dd(1,ix) - adt*work.dDu2dub(1:nz)./dzp;
    dd(2,ix) = dd(2,ix) - adt*(work.dDu2dua(1:nz) - work.dDu2dub(2:nzp))./dzp;
    dd(3,ix) = dd(3,ix) + adt*work.dDu2dua(2:nzp)./dzp;
    % Relabelling terms - only keep dependence on u2
    dd(2,ix) = dd(2,ix) - adt*M21.*relabel.duhat21du2 ...
                        + adt*M12.*relabel.duhat12du2;
    
    % Now solve the linear system
    % disp('*** zero u2 increments ***')
    inc_u2 = Ndiagsolveb(dd,rhsu);
    %inc_u2 = Ndiagsolvex(dd,rhsu);
    
    if conv_diag
        % For testing, compute predicted residuals using linearization
        % to compare with actual residuals ...
        yy = state_err.fluid(2).u;
        rr = - Ndiagmult(dd,yy);
        disp(' ')
        disp('Predicted Residual')
        disp(['res_u2      ',num2str(rr(krange))])
        zz = - adt*M21.*relabel.duhat21du1 + adt*M12.*relabel.duhat12du1;
        rrz = - zz.*state_err.fluid(1).u;
        rr = rr + rrz;
        disp(['with u1pert ',num2str(rr(krange))])
        %pause
    end
    
    % And increment u2
    state_new.fluid(2).u = state_new.fluid(2).u + inc_u2;
    if conv_diag
        disp(['inc_u2 = ',num2str(inc_u2(krange))])
    end
    
    % --------
    
    % Build linear system for v1 system
    ix = 1:nz;
    rhsu = res1v;
    
    dd = zeros(3,nz);
    % Tendency term
    dd(2,ix) = dd(2,ix) + m1;
    % Transport terms    
    dd(1,ix) = dd(1,ix) - adt*work.dFv1dvb(1:nz)./dzp;
    dd(2,ix) = dd(2,ix) - adt*(work.dFv1dva(1:nz) - work.dFv1dvb(2:nzp))./dzp;
    dd(3,ix) = dd(3,ix) + adt*work.dFv1dva(2:nzp)./dzp;
    % Diffusion terms
    dd(1,ix) = dd(1,ix) - adt*work.dDv1dvb(1:nz)./dzp;
    dd(2,ix) = dd(2,ix) - adt*(work.dDv1dva(1:nz) - work.dDv1dvb(2:nzp))./dzp;
    dd(3,ix) = dd(3,ix) + adt*work.dDv1dva(2:nzp)./dzp;
    % Relabelling terms - only keep dependence on v1
    dd(2,ix) = dd(2,ix) - adt*M12.*relabel.dvhat12dv1 ...
                        + adt*M21.*relabel.dvhat21dv1;
       
    % Now solve the linear system
    % disp('*** zero increments ***')
    inc_v1 = Ndiagsolveb(dd,rhsu);
    %inc_v1 = Ndiagsolvex(dd,rhsu);
    
    if conv_diag
        % For testing, compute predicted residuals using linearization
        % to compare with actual residuals ...
        yy = state_err.fluid(1).v;
        rr = - Ndiagmult(dd,yy);
        disp(' ')
        disp('Predicted Residual')
        disp(['res_v1    ' num2str(rr(krange))])
        zz = - adt*M12.*(relabel.dvhat12dv2) + adt*M21.*(relabel.dvhat21dv2);
        rrz = - zz.*state_err.fluid(2).v;
        rr = rr + rrz;
        disp(['with v2pert ',num2str(rr(krange))])
        %pause
    end
 
    % And increment v1
    state_new.fluid(1).v = state_new.fluid(1).v + inc_v1;
    if conv_diag
        disp(['inc_v1 = ',num2str(inc_v1(krange))])
    end
    
    % --------
    
    % Build linear system for v2 system
    ix = 1:nz;
    rhsu = res2v;
    
    dd = zeros(3,nz);
    % Tendency term
    if ~switches.c
        dd(2,ix) = dd(2,ix) + m2;
    else
        dd(2,ix) = dd(2,ix) + 0.01*m2;
    end
    % Transport terms    
    dd(1,ix) = dd(1,ix) - adt*work.dFv2dvb(1:nz)./dzp;
    dd(2,ix) = dd(2,ix) - adt*(work.dFv2dva(1:nz) - work.dFv2dvb(2:nzp))./dzp;
    dd(3,ix) = dd(3,ix) + adt*work.dFv2dva(2:nzp)./dzp;
    % Diffusion terms
    dd(1,ix) = dd(1,ix) - adt*work.dDv2dvb(1:nz)./dzp;
    dd(2,ix) = dd(2,ix) - adt*(work.dDv2dva(1:nz) - work.dDv2dvb(2:nzp))./dzp;
    dd(3,ix) = dd(3,ix) + adt*work.dDv2dva(2:nzp)./dzp;
    % Relabelling terms - only keep dependence on v2
    dd(2,ix) = dd(2,ix) - adt*M21.*relabel.dvhat21dv2 ...
                        + adt*M12.*relabel.dvhat12dv2;
          
    % Now solve the linear system
    % disp('*** zero v2 increments ***')
    inc_v2 = Ndiagsolveb(dd,rhsu);
    %inc_v2 = Ndiagsolvex(dd,rhsu);
    
    if conv_diag
        % For testing, compute predicted residuals using linearization
        % to compare with actual residuals ...
        yy = state_err.fluid(2).v;
        rr = - Ndiagmult(dd,yy);
        disp(' ')
        disp('Predicted Residual')
        disp(['res_v2      ',num2str(rr(krange))])
        zz = - adt*M21.*(relabel.dvhat21dv1) + adt*M12.*(relabel.dvhat12dv1);
        rrz = - zz.*state_err.fluid(1).v;
        rr = rr + rrz;
        disp(['with v1pert ',num2str(rr(krange))])
        %pause
    end
    
    % And increment v2
    %disp('*** No v2 inc ***')
    state_new.fluid(2).v = state_new.fluid(2).v + inc_v2;
    if conv_diag
        disp(['inc_v2 = ',num2str(inc_v2(krange))])
    end

    % --------
    
    % Linear system for all second moments
    build_2M_linear_system

    % Now build the right hand side
    ix = 1:12:12*nz-11;
    rhs2M(ix) = res1tke;
    ix = ix + 1;
    rhs2M(ix) = res2tke;
    ix = ix + 1;
    rhs2M(ix) = 0;
    ix = ix + 1;
    rhs2M(ix) = 0;
    ix = ix + 1;
    rhs2M(ix) = 0;
    ix = ix + 1;
    rhs2M(ix) = 0;
    ix = ix + 1;
    rhs2M(ix) = tend.fluid(1).mvareta.tot;
    ix = ix + 1;
    rhs2M(ix) = tend.fluid(2).mvareta.tot;
    ix = ix + 1;
    rhs2M(ix) = tend.fluid(1).mvarq.tot;
    ix = ix + 1;
    rhs2M(ix) = tend.fluid(2).mvarq.tot;
    ix = ix + 1;
    rhs2M(ix) = tend.fluid(1).mcovaretaq.tot;
    ix = ix + 1;
    rhs2M(ix) = tend.fluid(2).mcovaretaq.tot;
    
    % Now solve the linear system    
    xx = Ndiagsolveb(dd,rhs2M);
    %xx = Ndiagsolvex(dd,rhs2M);
    
    % and unpack the increments
    inc_tke1       = xx( 1:12:12*nz-11);
    inc_tke2       = xx( 2:12:12*nz-10);
    inc_weta1      = xx( 3:12:12*nz- 9);
    inc_weta2      = xx( 4:12:12*nz- 8);
    inc_wq1        = xx( 5:12:12*nz- 7);
    inc_wq2        = xx( 6:12:12*nz- 6);
    inc_vareta1    = xx( 7:12:12*nz- 5);
    inc_vareta2    = xx( 8:12:12*nz- 4);
    inc_varq1      = xx( 9:12:12*nz- 3);
    inc_varq2      = xx(10:12:12*nz- 2);
    inc_covaretaq1 = xx(11:12:12*nz- 1);
    inc_covaretaq2 = xx(12:12:12*nz);
    
%     if conv_diag
%         % For testing, compute predicted residuals using linearization
%         % to compare with actual residuals ...
%         % *** In order to do this we need tke and variance errors - see above ***
%         rr = - Ndiagmult(dd,zz);
%         % etc...
%     end
 
    % Increment tke1
    % disp('*** no tke1 inc ***')
    state_new.fluid(1).tke = state_new.fluid(1).tke + inc_tke1;
    inc_fix = max(constants.param.tke_min - state_new.fluid(1).tke,0);
    state_new.fluid(1).tke = state_new.fluid(1).tke + inc_fix;
    fixmtke1 = fixmtke1 + m1.*inc_fix;
    if conv_diag
        disp(['inc_tke1 = ',num2str(inc_tke1(krange))])
    end
    
    % Increment tke2
    % disp('*** no tke2 inc ***')
    state_new.fluid(2).tke = state_new.fluid(2).tke + inc_tke2;
    inc_fix = max(constants.param.tke_min - state_new.fluid(2).tke,0);
    state_new.fluid(2).tke = state_new.fluid(2).tke + inc_fix;
    fixmtke2 = fixmtke2 + m2.*inc_fix;
    if conv_diag
        disp(['inc_tke2 = ',num2str(inc_tke2(krange))])
        disp(['inc_fix  = ',num2str(inc_fix(krange))])
    end
       
    % Increment variances and covariances
% disp('*** bounded var decrements ***')
% disp('*** frozen covariances ***')
    state_new.fluid(1).vareta    = max(0.1*state_new.fluid(1).vareta,state_new.fluid(1).vareta + inc_vareta1);
    state_new.fluid(2).vareta    = max(0.1*state_new.fluid(2).vareta,state_new.fluid(2).vareta + inc_vareta2);
    state_new.fluid(1).varq      = max(0.1*state_new.fluid(1).varq,state_new.fluid(1).varq     + inc_varq1  );
    state_new.fluid(2).varq      = max(0.1*state_new.fluid(2).varq,state_new.fluid(2).varq     + inc_varq2  );
    state_new.fluid(1).covaretaq = state_new.fluid(1).covaretaq + inc_covaretaq1;
    state_new.fluid(2).covaretaq = state_new.fluid(2).covaretaq + inc_covaretaq2;
    
    % --------
    
    % Check for NaNs in increments
    % check_this = 'inc';
    % check_for_nans
    
    % Check for equality of increments when fluid 1 and 2
    % have identical properties
    % check_equal_increments
    
    % --------
        
    if min(state_new.fluid(1).m) < 0
        disp('m1 < 0')
        pause
    end
    if min(state_new.fluid(2).m) < 0
        disp('m2 < 0')
        pause
    end
    % pause
    
end
    
% Accumulate forcing and budget terms
accdt = adt;
accumulate
accumulate_fix

% diagnose_sg_flux