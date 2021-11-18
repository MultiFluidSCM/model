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
    wmean = sigma1w.*w1 + sigma2w.*w2;
    tke1_res = 3*0.5*(w1-wmean).^2;
    tke2_res = 3*0.5*(w2-wmean).^2;
    tke2_res = (abovep.*tke2_res(2:nzp) + belowp.*tke2_res(1:nz));
    velocity_scale = mix.tke1_factor*sqrt(tke1) + mix.tke2_factor*sqrt(tke2);
    % velocity_scale = mix.tke1_factor*sqrt(tke1) + mix.tke2_factor*sqrt(tke2_res);
    % velocity_scale = mix.tke1_factor*sqrt(tke1) + mix.tke2_factor*sqrt(0.5*tke2+0.5*tke2_res);
    r2 = min(velocity_scale./scales.L_plume, rdt);
    
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

% Default transfer coefficients
bdetrainw_mix = mix.bdetrainw;
bentrainw_mix = mix.bentrainw;
bdetrainq_mix = mix.bdetrainq;
bentrainq_mix = mix.bentrainq;
bdetraint_mix = mix.bdetraint;
bentraint_mix = mix.bentraint;
bdetrainw_mix_cloud = mix_cloud.bdetrainw;
bentrainw_mix_cloud = mix_cloud.bentrainw;
bdetrainq_mix_cloud = mix_cloud.bdetrainq;
bentrainq_mix_cloud = mix_cloud.bentrainq;
bdetraint_mix_cloud = mix_cloud.bdetraint;
bentraint_mix_cloud = mix_cloud.bentraint;

% Calculate the transfer (b) coefficients based on the transfer rate and the PDFs
if constants.param.mix.use_pdf
    deltaw   = w2-w1 + 1e-8*(abs(w2-w1) == 0);
    deltaq   = q2-q1 + 1e-8*(abs(q2-q1) == 0);
    deltaeta = eta2-eta1 + 1e-8*(abs(eta2-eta1) == 0);
    
    % Detrainment
    if mix.detrain
        w2std   = sqrt(tke2*2/3);
        q2std   = sqrt(state.fluid(2).varq);
        eta2std = sqrt(state.fluid(2).vareta);
        
        % Avoid division by zero
        w2std   = w2std   + 1e-8*(w2std   == 0);
        q2std   = q2std   + 1e-8*(q2std   == 0);
        eta2std = eta2std + 1e-8*(eta2std == 0);
        
        % Standard deviation at w-levels
        w2stdw   = weight_to_w(grid, w2std);
        q2stdw   = weight_to_w(grid, q2std);
        eta2stdw = weight_to_w(grid, eta2std);
        
        % Quantity which described how much of a pdf is transferred based on the transfer rate
        % See McIntyre 2020 Thesis, chapter 4.2
        pdf_factor12 = sqrt(2)*erfinv(2*dt*dM12dm2_mix - 1);
        
        % Part of PDF transferred is -infinity to cutoff below
        wcutoff12   = w2   + weight_to_w(grid, w2std   .* pdf_factor12);
        qcutoff12   = q2   + weight_to_w(grid, q2std   .* pdf_factor12);
        etacutoff12 = eta2 + weight_to_w(grid, eta2std .* pdf_factor12);
        
        % Mean of part of PDF transferred
        what12_mix = w2 - w2stdw .* exp(-0.5 * ((w2-wcutoff12)./w2stdw).^2) / sqrt(2*pi);
        qhat12_mix = q2 - q2stdw .* exp(-0.5 * ((q2-qcutoff12)./q2stdw).^2) / sqrt(2*pi);
        etahat12_mix = eta2 - eta2stdw .* exp(-0.5 * ((eta2-etacutoff12)./eta2stdw).^2) / sqrt(2*pi);
        
        bdetrainw_mix = (what12_mix - w1)./deltaw;
        bdetrainq_mix = (qhat12_mix - q1)./deltaq;
        bdetraint_mix = (etahat12_mix - eta1)./deltaeta;
        
        % Precaution to prevent too-extreme values
        bdetrainw_mix = min(2, max(-1, bdetrainw_mix));
        bdetrainq_mix = min(2, max(-1, bdetrainq_mix));
        bdetraint_mix = min(2, max(-1, bdetraint_mix));
        
        % Remove noise from the b-coefficients in places where no transfer is occuring.
        % If the transfer from dw/dz detrainment is very small, set the coefficient to 1.
        % If this is not done, the noise/spikes can be interpolated onto regions where transfer
        % is occuring which can cause the simulation to crash.
        filter = weight_to_w(grid, relabel.M12_mix) > 1e-4;
        % filter = weight_to_w(grid, relabel.M12_mix/(relabel.M12_mix + 1e-8*(relabel.M12_mix == 0)));
        bdetrainw_mix = bdetrainw_mix .* filter + (1-filter);
        bdetrainq_mix = bdetrainq_mix .* filter + (1-filter);
        bdetraint_mix = bdetraint_mix .* filter + (1-filter);
        
        bdetrainw_mix_cloud = bdetrainw_mix;
        bdetrainq_mix_cloud = bdetrainq_mix;
        bdetraint_mix_cloud = bdetraint_mix;
    end
    
    % Entrainment
    if mix.entrain
        w1std   = sqrt(tke1*2/3);
        q1std   = sqrt(state.fluid(1).varq);
        eta1std = sqrt(state.fluid(1).vareta);
        
        % Avoid division by zero
        w1std   = w1std   + 1e-8*(w1std   == 0);
        q1std   = q1std   + 1e-8*(q1std   == 0);
        eta1std = eta1std + 1e-8*(eta1std == 0);
        
        % Standard deviation at w-levels
        w1stdw   = weight_to_w(grid, w1std);
        q1stdw   = weight_to_w(grid, q1std);
        eta1stdw = weight_to_w(grid, eta1std);
        
        % Quantity which described how much of a pdf is transferred based on the transfer rate
        % See McIntyre 2020 Thesis, chapter 4.2
        pdf_factor21 = sqrt(2)*erfinv(2*dt*dM21dm1_mix - 1);
        
        % Part of PDF transferred is -infinity to cutoff below
        wcutoff21   = w1   - weight_to_w(grid, w1std   .* pdf_factor21);
        qcutoff21   = q1   - weight_to_w(grid, q1std   .* pdf_factor21);
        etacutoff21 = eta1 - weight_to_w(grid, eta1std .* pdf_factor21);
        
        % Mean of part of PDF transferred
        what21_mix = w1 + w1stdw .* exp(-0.5 * ((w1-wcutoff21)./w1stdw).^2) / sqrt(2*pi);
        qhat21_mix = q1 + q1stdw .* exp(-0.5 * ((q1-qcutoff21)./q1stdw).^2) / sqrt(2*pi);
        etahat21_mix = eta1 + eta1stdw .* exp(-0.5 * ((eta1-etacutoff21)./eta1stdw).^2) / sqrt(2*pi);
        
        bentrainw_mix = (what21_mix - w2)./-deltaw;
        bentrainq_mix = (qhat21_mix - q2)./-deltaq;
        bentraint_mix = (etahat21_mix - eta2)./-deltaeta;
        
        % Precaution to prevent too-extreme values
        bentrainw_mix = min(2, max(-1, bentrainw_mix));
        bentrainq_mix = min(2, max(-1, bentrainq_mix));
        bentraint_mix = min(2, max(-1, bentraint_mix));
        
        % Remove noise from the b-coefficients in places where no transfer is occuring.
        % If the transfer from dw/dz detrainment is very small, set the coefficient to 1.
        % If this is not done, the noise/spikes can be interpolated onto regions where transfer
        % is occuring which can cause the simulation to crash.
        filter = weight_to_w(grid, relabel.M21_mix) > 1e-4;
        % filter = weight_to_w(grid, relabel.M12_mix/(relabel.M12_mix + 1e-8*(relabel.M12_mix == 0)));
        bentrainw_mix = bentrainw_mix .* filter + (1-filter);
        bentrainq_mix = bentrainq_mix .* filter + (1-filter);
        bentraint_mix = bentraint_mix .* filter + (1-filter);
        
        bentrainw_mix_cloud = bentrainw_mix;
        bentrainq_mix_cloud = bentrainq_mix;
        bentraint_mix_cloud = bentraint_mix;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transfer coefficients in the boundary layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Entrained and detrained values of eta
[relabel.etahat12_mix_dry,relabel.detahat12deta1_mix_dry,relabel.detahat12deta2_mix_dry,...
 relabel.etahat21_mix_dry,relabel.detahat21deta1_mix_dry,relabel.detahat21deta2_mix_dry] ...
    = findqhat(eta1, eta2, bentraint_mix, bdetraint_mix);

% Entrained and detrained values of q
[relabel.qhat12_mix_dry,relabel.dqhat12dq1_mix_dry,relabel.dqhat12dq2_mix_dry,...
 relabel.qhat21_mix_dry,relabel.dqhat21dq1_mix_dry,relabel.dqhat21dq2_mix_dry] ...
    = findqhat(q1, q2, bentrainq_mix, bdetrainq_mix);

% Entrained and detrained values of velocities u, v and w
[relabel.uhat12_mix_dry,relabel.duhat12du1_mix_dry,relabel.duhat12du2_mix_dry,...
 relabel.uhat21_mix_dry,relabel.duhat21du1_mix_dry,relabel.duhat21du2_mix_dry] ...
    = findqhat(u1, u2, mix.bentrainu, mix.bdetrainu);

[relabel.vhat12_mix_dry,relabel.dvhat12dv1_mix_dry,relabel.dvhat12dv2_mix_dry,...
 relabel.vhat21_mix_dry,relabel.dvhat21dv1_mix_dry,relabel.dvhat21dv2_mix_dry] ...
    = findqhat(v1, v2, mix.bentrainu, mix.bdetrainu);

[relabel.what12_mix_dry,relabel.dwhat12dw1_mix_dry,relabel.dwhat12dw2_mix_dry,...
 relabel.what21_mix_dry,relabel.dwhat21dw1_mix_dry,relabel.dwhat21dw2_mix_dry] ...
    = findqhat(w1, w2, bentrainw_mix, bdetrainw_mix);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transfer coefficients in the cloud layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Entrained and detrained values of eta
[relabel.etahat12_mix_cloud,relabel.detahat12deta1_mix_cloud,relabel.detahat12deta2_mix_cloud,...
 relabel.etahat21_mix_cloud,relabel.detahat21deta1_mix_cloud,relabel.detahat21deta2_mix_cloud] ...
    = findqhat(eta1, eta2, bentraint_mix_cloud, bdetraint_mix_cloud);

% Entrained and detrained values of q
[relabel.qhat12_mix_cloud,relabel.dqhat12dq1_mix_cloud,relabel.dqhat12dq2_mix_cloud,...
 relabel.qhat21_mix_cloud,relabel.dqhat21dq1_mix_cloud,relabel.dqhat21dq2_mix_cloud] ...
    = findqhat(q1, q2, bentrainq_mix_cloud, bdetrainq_mix_cloud);

% Entrained and detrained values of velocities u, v and w
[relabel.uhat12_mix_cloud,relabel.duhat12du1_mix_cloud,relabel.duhat12du2_mix_cloud,...
 relabel.uhat21_mix_cloud,relabel.duhat21du1_mix_cloud,relabel.duhat21du2_mix_cloud] ...
    = findqhat(u1, u2, mix_cloud.bentrainu, mix_cloud.bdetrainu);

[relabel.vhat12_mix_cloud,relabel.dvhat12dv1_mix_cloud,relabel.dvhat12dv2_mix_cloud,...
 relabel.vhat21_mix_cloud,relabel.dvhat21dv1_mix_cloud,relabel.dvhat21dv2_mix_cloud] ...
    = findqhat(v1, v2, mix_cloud.bentrainu, mix_cloud.bdetrainu);

[relabel.what12_mix_cloud,relabel.dwhat12dw1_mix_cloud,relabel.dwhat12dw2_mix_cloud,...
 relabel.what21_mix_cloud,relabel.dwhat21dw1_mix_cloud,relabel.dwhat21dw2_mix_cloud] ...
    = findqhat(w1, w2, bentrainw_mix_cloud, bdetrainw_mix_cloud);


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