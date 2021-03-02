% Instability source entrainment where stratification is unstable
% *** Put tuneable coefficients in constants.params ***
relabel.M12_instab = instab.detrain * 0.  * m2 .* n2pos;
relabel.M21_instab = instab.entrain * 0.2 * m1 .* n1neg;

dM12dm2_instab = instab.detrain * 0.  * n2pos;
dM21dm1_instab = instab.entrain * 0.2 * n1neg;

% Entrained and detrained values of eta
[relabel.etahat12_instab,relabel.detahat12deta1_instab,relabel.detahat12deta2_instab,...
 relabel.etahat21_instab,relabel.detahat21deta1_instab,relabel.detahat21deta2_instab] ...
    = findqhat(eta1, eta2, instab.bentraint, instab.bdetraint);

% Entrained and detrained values of q
[relabel.qhat12_instab,relabel.dqhat12dq1_instab,relabel.dqhat12dq2_instab,...
 relabel.qhat21_instab,relabel.dqhat21dq1_instab,relabel.dqhat21dq2_instab] ...
    = findqhat(q1, q2, instab.bentrainq, instab.bdetrainq);

% Entrained and detrained values of velocities u, v and w
[relabel.uhat12_instab,relabel.duhat12du1_instab,relabel.duhat12du2_instab,...
 relabel.uhat21_instab,relabel.duhat21du1_instab,relabel.duhat21du2_instab] ...
    = findqhat(u1, u2, instab.bentrainu, instab.bdetrainu);

[relabel.vhat12_instab,relabel.dvhat12dv1_instab,relabel.dvhat12dv2_instab,...
 relabel.vhat21_instab,relabel.dvhat21dv1_instab,relabel.dvhat21dv2_instab] ...
    = findqhat(v1, v2, instab.bentrainu, instab.bdetrainu);

[relabel.what12_instab,relabel.dwhat12dw1_instab,relabel.dwhat12dw2_instab,...
 relabel.what21_instab,relabel.dwhat21dw1_instab,relabel.dwhat21dw2_instab] ...
    = findqhat(w1, w2, instab.bentrainw, instab.bdetrainw);

% relabel.what12_instab = min(relabel.what21_instab, 0);
% relabel.what21_instab = max(relabel.what21_instab, 0);