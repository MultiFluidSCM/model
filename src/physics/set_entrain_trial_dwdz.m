% Dynamical entrainment/detrainment based on vertical velocity convergence
rate_sort12 = 0*m2;
rate_sort21 = 0*m1;

if ischeme == 4
    % Maximum area fraction for downdraft (w < 0) regions from LES
    sigmad = 0.5;
    % Fraction of fluid 1 that is a downdraft
    frac1_down = min(sigma1, sigmad);
    % Fraction of transfer that actually enters fluid 2
    frac1_transfer = 2*min(sigma2, 0.5);
    
    rate_sort12 = min(max(0, -dw2dz), rdt);
    rate_sort21 = min(max(0, -dw1dz), rdt) .* frac1_down .* frac1_transfer;
end

relabel.M12_dwdz = dwdz.detrain * m2.*rate_sort12;
relabel.M21_dwdz = dwdz.entrain * m1.*rate_sort21;

dM12dm1_dwdz = 0 * dwdz.detrain * rate_sort12;
dM12dm2_dwdz =     dwdz.detrain * rate_sort12;

dM21dm1_dwdz =     dwdz.entrain * rate_sort21;
dM21dm2_dwdz = 0 * dwdz.entrain * rate_sort21;

% Prevent double counting from instab entrainment 
% relabel.M21_dwdz = max(relabel.M21_dwdz-relabel.M21_instab, 0);
% dM21dm1_dwdz = max(dM21dm1_dwdz-dM21dm1_instab, 0);

% Entrained and detrained values of eta
[relabel.etahat12_dwdz,relabel.detahat12deta1_dwdz,relabel.detahat12deta2_dwdz,...
 relabel.etahat21_dwdz,relabel.detahat21deta1_dwdz,relabel.detahat21deta2_dwdz] ...
    = findqhat(eta1, eta2, dwdz.bentraint, dwdz.bdetraint);

% Entrained and detrained values of q
[relabel.qhat12_dwdz,relabel.dqhat12dq1_dwdz,relabel.dqhat12dq2_dwdz,...
 relabel.qhat21_dwdz,relabel.dqhat21dq1_dwdz,relabel.dqhat21dq2_dwdz] ...
    = findqhat(q1, q2, dwdz.bentrainq, dwdz.bdetrainq);

% Entrained and detrained values of velocities u, v and w
[relabel.uhat12_dwdz,relabel.duhat12du1_dwdz,relabel.duhat12du2_dwdz,...
 relabel.uhat21_dwdz,relabel.duhat21du1_dwdz,relabel.duhat21du2_dwdz] ...
    = findqhat(u1, u2, dwdz.bentrainu, dwdz.bdetrainu);

[relabel.vhat12_dwdz,relabel.dvhat12dv1_dwdz,relabel.dvhat12dv2_dwdz,...
 relabel.vhat21_dwdz,relabel.dvhat21dv1_dwdz,relabel.dvhat21dv2_dwdz] ...
    = findqhat(v1, v2, dwdz.bentrainu, dwdz.bdetrainu);

[relabel.what12_dwdz,relabel.dwhat12dw1_dwdz,relabel.dwhat12dw2_dwdz,...
 relabel.what21_dwdz,relabel.dwhat21dw1_dwdz,relabel.dwhat21dw2_dwdz] ...
    = findqhat(w1, w2, dwdz.bentrainw, dwdz.bdetrainw);
    
relabel.what12_dwdz = 0*relabel.what12_dwdz;
relabel.dwhat12dw1_dwdz = 0*relabel.dwhat12dw1_dwdz;
relabel.dwhat12dw2_dwdz = 0*relabel.dwhat12dw2_dwdz;

relabel.what21_dwdz = 0*relabel.what21_dwdz;
relabel.dwhat21dw1_dwdz = 0*relabel.dwhat21dw1_dwdz;
relabel.dwhat21dw2_dwdz = 0*relabel.dwhat21dw2_dwdz;


% relabel.what12_dwdz = min(relabel.what21_dwdz, 0);
% relabel.what21_dwdz = max(relabel.what21_dwdz, 0);