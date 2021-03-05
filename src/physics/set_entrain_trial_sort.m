% APDF-based formulation for detrainment, sorting pdf
relabel.w2bar = w2bar;
relabel.wstd = wstd;

% 1/sqrt(2) * normalized variable at cut of pdf on p-levels
chi_cut = - rroot2*w2bar./wstd;

% and w-levels (bounded to prevent divide by zero below)
% chi_cut_w = max(-rroot2*w2./wstdw,-4);
chi_cut_w = -rroot2*w2./wstdw;

% Fraction of pdf cut off i.e. with w < 0
relabel.frac = 0.5*(1 + erf(chi_cut));

% Rate at which it is cut off
% rate_sort = 5*(wstd./scales.L_plume);
if ischeme == 1
    rate_sort = min(10*max(0,-dw2dz),rdt);
elseif ischeme == 3 | ischeme == 4
    % rate_sort = min(20*max(0,-dw2dz),rdt);
    % rate_sort = sign(n1pos).*min(20*rate_mix,rdt);  % Strictly this is scheme 7 or 8
    % rate_sort = (buoybar < 0).*min(20*rate_mix,rdt);  % 
    % ss3 = 0.5*(1 - tanh( 2*buoy./(bstd + eps) + 0.4));  % Smooth switch 10b
    % ss3 = 0.5*(1 - tanh( 2*buoy./(bstd + eps) + 0.8));  % Smooth switch 10d
    % ss3 = 0.5*(1 - tanh( 5*buoy./(bstd + eps) + 2));  % Smooth switch 10c
    bnorm = weight_to_w(grid,tke2./scales.L_turb2);
    ss3 = 0.5*(1 - tanh( 3*buoy./bnorm + 1.2));  % Smooth switch 10e
    relabel.ss3bar = abovep.*ss3(2:nzp) + belowp.*ss3(1:nz);
    rate_sort = relabel.ss3bar.*min(10*rate_mix,rdt);  % 
    % rate_buoy = max(-buoybar,0)./wstd;
    % rate_sort = min(20*rate_buoy,rdt);
    
    % rate_sort = rate_sort + max(0,-dw2dz);
    % rate_sort = min(rate_sort,rdt);
    
    % [temp, cbase_index] = min(abs(grid.zw-2200));
    % for k = cbase_index:length(rate_sort)
        % rate_sort(k) = 0;
    % end
else
    rate_sort = 0;
end
% rate_sort = 0;

% Resulting detrainment rate
relabel.M12_sort = sort.detrain * rate_sort.*m2.*relabel.frac;
relabel.M21_sort = sort.entrain * 0 * m1;


% w, eta and q of detrained fluid
% chi_hat = - sqrt(2/pi)*exp(-chi_cut_w.^2)./(1 + erf(chi_cut_w));
for k = 1:nzp
    [chi_hat(k),deriv(k)] = compute_chi_hat(chi_cut_w(k));
end

% Detrained fluid properties
relabel.uhat12_sort   = u2;
relabel.vhat12_sort   = v2;
relabel.what12_sort   = w2   + wstdw .*chi_hat;
relabel.etahat12_sort = eta2 + etastd.*chi_hat;
relabel.qhat12_sort   = q2   + qstd  .*chi_hat;

% Entrained fluid properties
relabel.uhat21_sort   = u1;
relabel.vhat21_sort   = v1;
relabel.what21_sort   = w1;
relabel.etahat21_sort = eta1;
relabel.qhat21_sort   = q1;

% Derivative of what12_sort wrt wstdw
relabel.what_deriv = deriv;


% Property derivatives. Assume transferred properties close to fluid mean values.
relabel.detahat12deta1_sort = zeros(1,nzp);
relabel.detahat12deta2_sort = ones(1,nzp);
relabel.detahat21deta1_sort = ones(1,nzp);
relabel.detahat21deta2_sort = zeros(1,nzp);

relabel.dqhat12dq1_sort = zeros(1,nzp);
relabel.dqhat12dq2_sort = ones(1,nzp);
relabel.dqhat21dq1_sort = ones(1,nzp);
relabel.dqhat21dq2_sort = zeros(1,nzp);

relabel.duhat12du1_sort = zeros(1,nz);
relabel.duhat12du2_sort = ones(1,nz);
relabel.duhat21du1_sort = ones(1,nz);
relabel.duhat21du2_sort = zeros(1,nz);

relabel.dvhat12dv1_sort = zeros(1,nz);
relabel.dvhat12dv2_sort = ones(1,nz);
relabel.dvhat21dv1_sort = ones(1,nz);
relabel.dvhat21dv2_sort = zeros(1,nz);

relabel.dwhat12dw1_sort = zeros(1,nzp);
relabel.dwhat12dw2_sort = ones(1,nzp);
relabel.dwhat21dw1_sort = ones(1,nzp);
relabel.dwhat21dw2_sort = zeros(1,nzp);