% Relaxation to homogeneous state with background sigma
% Ratio of sigma2 to reference value
rr = sigma2/sigma20;
% Merging factor to go smoothly between fast relaxation rate r1 for rr < 1
% to turbulent rate r2 for rr > 2
rr = min(2, max(1,rr) );
qq = (rr - 2).*(rr - 2).*(1 - (rr - 1).*(rr - 1));
r1 = rdt;
if ischeme == 0
    r2 = n1pos;
elseif ischeme == 1
    r2 = min(0.5*sqrt(tke2)./scales.L_plume,rdt);
elseif ischeme == 3 | ischeme == 4
    % r2 = min(0.25*sqrt(tke2)./scales.L_plume,rdt);
    r2 = min(mix.tke_factor*sqrt(tke1)./scales.L_plume,rdt);
    % r2 = min(0.25*sqrt(tke2)./grid.zp,rdt);
    % r2 = min(0.25*sqrt(sigma2.*tke2)./scales.L_plume,rdt);
    % r2 = min(0.25*sqrt(tke1)./scales.L_plume,rdt);
    % r2 = min(0.25*sqrt(tke2)./scales.L_turb2,rdt);
    % r2 = min(0.25*sqrt(tke2)./(sigma1.*scales.L_turb1+sigma2.*scales.L_turb2),rdt);
    % r2 = min(0.25*sqrt(tke1+tke2)./(scales.L_turb1+scales.L_turb2),rdt);
    
    % Caltech mixing formulation
    % zFluid2Max = 100;
    % for i=1:length(sigma2)
        % if (sigma2(i)-sigma20) > 0
            % zFluid2Max = max(100, grid.zp(i));
        % end
    % end
    % r2 = min(0.15*sqrt(tke1)/zFluid2Max,rdt);
    % r2 = min(sqrt(tke1)/zFluid2Max,rdt);
    % r2 = min(mix.tke_factor*sqrt(tke2)/zFluid2Max,rdt);
    % r2 = min(sqrt(sigma1.*tke1+sigma2.*tke2)/zFluid2Max,rdt);
else
    disp('unknown scheme in set_entrain_trial')
    pause
end
rate_mix = r2 + (r1 - r2).*qq;
relabel.rate_mix = rate_mix; % Just for diagnostics
% Reference profiles of  m1 and m2
m2ref = sigma20*rho;
m1ref = rho - m2ref;
% Factors needed to compute mixing
m1fac = 2*mf1 - sigma10;
m2fac = 2*mf2 - sigma20;
% There are three cases:
% m1fac < 0, 0 < m1fac < 1, 1 < m1fac
case1 = m1fac < 0;
case3 = 1 < m1fac;
case2 = 1 - (case1 + case3);

% Bound to prevent negative relabelling rates
m1fac = min(1, max(0,m1fac) );
m2fac = min(1, max(0,m2fac) );

% Relabelling rates
relabel.M21_mix = mix.entrain * rate_mix.*m2.*m1fac;
relabel.M12_mix = mix.detrain * rate_mix.*m1.*m2fac;

% [temp, cbase_index] = min(abs(grid.zw-2000));
% for k = cbase_index:length(relabel.M12_mix)
    % relabel.M12_mix(k) = relabel.M12_mix(k)*2000/grid.zw(k);
% end

% Derivatives *** Need to think carefully about these ***
dM21dm1_mix = mix.entrain * (case1.*zeros(1,nz) + case2.*2.*rate_mix.*mf2.*mf2            + case3.*zeros(1,nz));
dM21dm2_mix = mix.entrain * (case1.*zeros(1,nz) + case2.*rate_mix.*(2*mf1.*mf1 - sigma10) + case3.*rate_mix);
dM12dm1_mix = mix.detrain * (case1.*rate_mix    + case2.*rate_mix.*(2*mf2.*mf2 - sigma20) + case3.*zeros(1,nz));
dM12dm2_mix = mix.detrain * (case1.*zeros(1,nz) + case2.*2.*rate_mix.*mf1.*mf1            + case3.*zeros(1,nz));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transfer coefficients in the boundary layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Entrained and detrained values of eta
[relabel.etahat12_mix_dry,relabel.detahat12deta1_mix_dry,relabel.detahat12deta2_mix_dry,...
 relabel.etahat21_mix_dry,relabel.detahat21deta1_mix_dry,relabel.detahat21deta2_mix_dry] ...
    = findqhat(eta1, eta2, mix.bentraint, mix.bdetraint);

% Entrained and detrained values of q
[relabel.qhat12_mix_dry,relabel.dqhat12dq1_mix_dry,relabel.dqhat12dq2_mix_dry,...
 relabel.qhat21_mix_dry,relabel.dqhat21dq1_mix_dry,relabel.dqhat21dq2_mix_dry] ...
    = findqhat(q1, q2, mix.bentrainq, mix.bdetrainq);

% Entrained and detrained values of velocities u, v and w
[relabel.uhat12_mix_dry,relabel.duhat12du1_mix_dry,relabel.duhat12du2_mix_dry,...
 relabel.uhat21_mix_dry,relabel.duhat21du1_mix_dry,relabel.duhat21du2_mix_dry] ...
    = findqhat(u1, u2, mix.bentrainu, mix.bdetrainu);

[relabel.vhat12_mix_dry,relabel.dvhat12dv1_mix_dry,relabel.dvhat12dv2_mix_dry,...
 relabel.vhat21_mix_dry,relabel.dvhat21dv1_mix_dry,relabel.dvhat21dv2_mix_dry] ...
    = findqhat(v1, v2, mix.bentrainu, mix.bdetrainu);

[relabel.what12_mix_dry,relabel.dwhat12dw1_mix_dry,relabel.dwhat12dw2_mix_dry,...
 relabel.what21_mix_dry,relabel.dwhat21dw1_mix_dry,relabel.dwhat21dw2_mix_dry] ...
    = findqhat(w1, w2, mix.bentrainw, mix.bdetrainw);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transfer coefficients in the cloud layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Entrained and detrained values of eta
[relabel.etahat12_mix_cloud,relabel.detahat12deta1_mix_cloud,relabel.detahat12deta2_mix_cloud,...
 relabel.etahat21_mix_cloud,relabel.detahat21deta1_mix_cloud,relabel.detahat21deta2_mix_cloud] ...
    = findqhat(eta1, eta2, mix_cloud.bentraint, mix_cloud.bdetraint);

% Entrained and detrained values of q
[relabel.qhat12_mix_cloud,relabel.dqhat12dq1_mix_cloud,relabel.dqhat12dq2_mix_cloud,...
 relabel.qhat21_mix_cloud,relabel.dqhat21dq1_mix_cloud,relabel.dqhat21dq2_mix_cloud] ...
    = findqhat(q1, q2, mix_cloud.bentrainq, mix_cloud.bdetrainq);

% Entrained and detrained values of velocities u, v and w
[relabel.uhat12_mix_cloud,relabel.duhat12du1_mix_cloud,relabel.duhat12du2_mix_cloud,...
 relabel.uhat21_mix_cloud,relabel.duhat21du1_mix_cloud,relabel.duhat21du2_mix_cloud] ...
    = findqhat(u1, u2, mix_cloud.bentrainu, mix_cloud.bdetrainu);

[relabel.vhat12_mix_cloud,relabel.dvhat12dv1_mix_cloud,relabel.dvhat12dv2_mix_cloud,...
 relabel.vhat21_mix_cloud,relabel.dvhat21dv1_mix_cloud,relabel.dvhat21dv2_mix_cloud] ...
    = findqhat(v1, v2, mix_cloud.bentrainu, mix_cloud.bdetrainu);

[relabel.what12_mix_cloud,relabel.dwhat12dw1_mix_cloud,relabel.dwhat12dw2_mix_cloud,...
 relabel.what21_mix_cloud,relabel.dwhat21dw1_mix_cloud,relabel.dwhat21dw2_mix_cloud] ...
    = findqhat(w1, w2, mix_cloud.bentrainw, mix_cloud.bdetrainw);


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
    % [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(pbar,Tw1(k),q1(k),constants.therm);
    % vapour1(k) = (1 - q1(k))*(1 - a)/a;
    % liquid1(k) = q1(k) - vapour1(k);
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(pbar,Tw2(k),q2(k),constants.therm);
    vapour2(k) = (1 - q2(k))*(1 - a)/a;
    liquid2(k) = q2(k) - vapour2(k);
end

% Smooth switch for weather a cloud is present or not
cloud = 0.5*(1 + tanh( 2e3*(liquid2-1e-5) ));
cloudp = grid.abovep.*cloud(2:nzp) + grid.belowp.*cloud(1:nz);

% Combine transfer coefficients
relabel.etahat12_mix       = (1-cloud).*relabel.etahat12_mix_dry       + cloud.*relabel.etahat12_mix_cloud;
relabel.etahat21_mix       = (1-cloud).*relabel.etahat21_mix_dry       + cloud.*relabel.etahat21_mix_cloud;
relabel.detahat12deta1_mix = (1-cloud).*relabel.detahat12deta1_mix_dry + cloud.*relabel.detahat12deta1_mix_cloud;
relabel.detahat21deta1_mix = (1-cloud).*relabel.detahat21deta1_mix_dry + cloud.*relabel.detahat21deta1_mix_cloud;
relabel.detahat12deta2_mix = (1-cloud).*relabel.detahat12deta2_mix_dry + cloud.*relabel.detahat12deta2_mix_cloud;
relabel.detahat21deta2_mix = (1-cloud).*relabel.detahat21deta2_mix_dry + cloud.*relabel.detahat21deta2_mix_cloud;

relabel.qhat12_mix     = (1-cloud).*relabel.qhat12_mix_dry     + cloud.*relabel.qhat12_mix_cloud;
relabel.qhat21_mix     = (1-cloud).*relabel.qhat21_mix_dry     + cloud.*relabel.qhat21_mix_cloud;
relabel.dqhat12dq1_mix = (1-cloud).*relabel.dqhat12dq1_mix_dry + cloud.*relabel.dqhat12dq1_mix_cloud;
relabel.dqhat21dq1_mix = (1-cloud).*relabel.dqhat21dq1_mix_dry + cloud.*relabel.dqhat21dq1_mix_cloud;
relabel.dqhat12dq2_mix = (1-cloud).*relabel.dqhat12dq2_mix_dry + cloud.*relabel.dqhat12dq2_mix_cloud;
relabel.dqhat21dq2_mix = (1-cloud).*relabel.dqhat21dq2_mix_dry + cloud.*relabel.dqhat21dq2_mix_cloud;

relabel.uhat12_mix     = (1-cloudp).*relabel.uhat12_mix_dry     + cloudp.*relabel.uhat12_mix_cloud;
relabel.uhat21_mix     = (1-cloudp).*relabel.uhat21_mix_dry     + cloudp.*relabel.uhat21_mix_cloud;
relabel.duhat12du1_mix = (1-cloudp).*relabel.duhat12du1_mix_dry + cloudp.*relabel.duhat12du1_mix_cloud;
relabel.duhat21du1_mix = (1-cloudp).*relabel.duhat21du1_mix_dry + cloudp.*relabel.duhat21du1_mix_cloud;
relabel.duhat12du2_mix = (1-cloudp).*relabel.duhat12du2_mix_dry + cloudp.*relabel.duhat12du2_mix_cloud;
relabel.duhat21du2_mix = (1-cloudp).*relabel.duhat21du2_mix_dry + cloudp.*relabel.duhat21du2_mix_cloud;

relabel.vhat12_mix     = (1-cloudp).*relabel.vhat12_mix_dry     + cloudp.*relabel.vhat12_mix_cloud;
relabel.vhat21_mix     = (1-cloudp).*relabel.vhat21_mix_dry     + cloudp.*relabel.vhat21_mix_cloud;
relabel.dvhat12dv1_mix = (1-cloudp).*relabel.dvhat12dv1_mix_dry + cloudp.*relabel.dvhat12dv1_mix_cloud;
relabel.dvhat21dv1_mix = (1-cloudp).*relabel.dvhat21dv1_mix_dry + cloudp.*relabel.dvhat21dv1_mix_cloud;
relabel.dvhat12dv2_mix = (1-cloudp).*relabel.dvhat12dv2_mix_dry + cloudp.*relabel.dvhat12dv2_mix_cloud;
relabel.dvhat21dv2_mix = (1-cloudp).*relabel.dvhat21dv2_mix_dry + cloudp.*relabel.dvhat21dv2_mix_cloud;

relabel.what12_mix     = (1-cloud).*relabel.what12_mix_dry     + cloud.*relabel.what12_mix_cloud;
relabel.what21_mix     = (1-cloud).*relabel.what21_mix_dry     + cloud.*relabel.what21_mix_cloud;
relabel.dwhat12dw1_mix = (1-cloud).*relabel.dwhat12dw1_mix_dry + cloud.*relabel.dwhat12dw1_mix_cloud;
relabel.dwhat21dw1_mix = (1-cloud).*relabel.dwhat21dw1_mix_dry + cloud.*relabel.dwhat21dw1_mix_cloud;
relabel.dwhat12dw2_mix = (1-cloud).*relabel.dwhat12dw2_mix_dry + cloud.*relabel.dwhat12dw2_mix_cloud;
relabel.dwhat21dw2_mix = (1-cloud).*relabel.dwhat21dw2_mix_dry + cloud.*relabel.dwhat21dw2_mix_cloud;


% Force transferred velocities to have closer values to their new fluid
% relabel.what12_mix = min(relabel.what21_mix, 0);
% relabel.what21_mix = max(relabel.what21_mix, 0);