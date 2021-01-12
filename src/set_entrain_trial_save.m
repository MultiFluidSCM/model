function relabel = set_entrain_trial(grid, state, buoy, eos, scales, kdiffw2, constants, dt)

% Set mass entrainment and detrainment rates
% M_ij = ( m_i m_j / rho_i ) * rate_ij
%
% the values of all entrained fields what, etahat, qhat, uhat, vhat
%
% and their derivatives

% First unpack the fields needed
nz = grid.nz;
nzp = nz + 1;
zp = grid.zp;
abovep = grid.abovep;
belowp = grid.belowp;
m1   = state.fluid(1).m;
m2   = state.fluid(2).m;
eta1 = state.fluid(1).eta;
eta2 = state.fluid(2).eta;
w1   = state.fluid(1).w;
w2   = state.fluid(2).w;
q1   = state.fluid(1).q;
q2   = state.fluid(2).q;
u1   = state.fluid(1).u;
u2   = state.fluid(2).u;
v1   = state.fluid(1).v;
v2   = state.fluid(2).v;
tke2 = state.fluid(2).tke;
bentraint = constants.param.bentraint;
bentrainq = constants.param.bentrainq;
bentrainw = constants.param.bentrainw;
bentrainu = constants.param.bentrainu;
bdetraint = constants.param.bdetraint;
bdetrainq = constants.param.bdetrainq;
bdetrainw = constants.param.bdetrainw;
bdetrainu = constants.param.bdetrainu;

sigma1 = eos.sigma1;
sigma2 = eos.sigma2;

% Background sigma2 when nothing is happening
sigma00 = constants.param.sigma00;

% To bound relaxation rate
rdt = 1/dt;

% zstar = scales.zstar;
% wstar = scales.wstar;

% Buoyancy on p levels
buoybar = (abovep.*buoy(2:nzp) + belowp.*buoy(1:nz));

% Mean density
rho = m1 + m2;

% Mass fractions
mf1 = m1./rho;
mf2 = m2./rho;


% % Ideal mass fraction
% cflim = 0.01;
%disp('** reduced E D above BL **')
%th = tanh((zp - scales.zstar)./(0.2*scales.zstar));
%factor = 0.55 - 0.45*th;
% ideal = 0.5*(1 + th)*(constants.param.confrac - cflim) + cflim;
% ideal = ideal.*((zp + constants.param.zrough)/zstar).^(-1/6);

% % Check for descent
% %wdown = min(w2(1:nz),w2(2:nzp));
% w2bar = abovep.*w2(2:nzp) + belowp.*w2(1:nz);
    
% % Rate
% %bgrate    =  0.1*max(wstar/zstar,0.0001);
% bgrate    =  0.2*0.12*max(wstar/zstar,0.0001);
% ratescale = 30.0*max(wstar/zstar,0.0001);


% Now set entrainment and detrainment rates
% and their derivatives for linearization
for k = 1:nz

    % *** Do we need a different buoyancy threshold here?
    if eos.nsq1(k) > 0 % | sigma2(k) < sigma00 % & buoybar(k) <= 0
        
        % Relax to background sigma and homogenized state
        
        % Time scale - consider using a fraction of rdt

        % Factor to make a smooth transition between fast relaxation
        % when sigma2 < sigma00 and relaxation at rate sqrt(-nsq1)
        % when sigma2 > 2*sigma00
        rr = sigma2(k)/sigma00;
        r1 = rdt;
        r2 = min(sqrt(eos.nsq1(k)),rdt);
        if rr < 1
            rate = r1;
        elseif rr > 2
            rate = r2;
        else
            qq = (rr - 2)*(rr - 2)*(1 - (rr - 1)*(rr - 1));
            rate = r2 + (r1 - r2)*qq;
        end
        %rate = rate*factor(k);
        
        % Reference m1 and m2
        m2ref = sigma00*rho(k);
        m1ref = rho(k) - m2ref;
        % Avoid negative transfers
        if m1(k) < 0.5*m1ref
            M21(k) = 0;
            M12(k) = m1(k)*rate;
            dM21dm1(k) = 0;
            dM21dm2(k) = 0;
            dM12dm1(k) = rate;
            dM12dm2(k) = 0;
        elseif m2(k) < 0.5*m2ref
            M12(k) = 0;
            M21(k) = m2(k)*rate;
            dM21dm1(k) = 0;
            dM21dm2(k) = rate;
            dM12dm1(k) = 0;
            dM12dm2(k) = 0;
        else
            M12(k) = sigma1(k)*(2*m2(k) - m2ref)*rate;
            M21(k) = sigma2(k)*(2*m1(k) - m1ref)*rate;
            dM21dm1(k) = 2*sigma2(k)*rate;
            dM21dm2(k) = (2*m1(k) - m1ref)*rate/rho(k);
            dM12dm1(k) = (2*m2(k) - m2ref)*rate/rho(k);
            dM12dm2(k) = 2*sigma1(k)*rate;
        end
        
    else
        
        % Source entrainment in unstable stratification
        % *** Put tuneable coefficients in constants.params ***
        rate = 0.2*sqrt(-eos.nsq1(k));
        M12(k) = 0;
        M21(k) = m1(k)*rate;
        dM21dm1(k) = rate;
        dM21dm2(k) = 0;
        dM12dm1(k) = 0;
        dM12dm2(k) = 0;
        
    end
    
    % Other derivatives
    dM21dw1(k) = 0;
    dM21dw2(k) = 0;
    dM21deta1(k) = 0;
    dM21deta2(k) = 0;
    dM21dq1(k) = 0;
    dM21dq2(k) = 0;

    dM12dw1(k) = 0;
    dM12dw2(k) = 0;
    dM12deta1(k) = 0;
    dM12deta2(k) = 0;
    dM12dq1(k) = 0;
    dM12dq2(k) = 0;
    
    % Trial scheme:
    % Relax to background sigma and homogenized state
    % on turbulent time scale

    % Factor to make a smooth transition between fast relaxation
    % when sigma2 < sigma00 and relaxation at turbulent rate
    % when sigma2 > 2*sigma00
    rr = sigma2(k)/sigma00;
    r1 = rdt;
    r2 = min(0.5*sqrt(tke2(k))/scales.L_plume(k),rdt);
    if rr < 1
        rate = r1;
    elseif rr > 2
        rate = r2;
    else
        qq = (rr - 2)*(rr - 2)*(1 - (rr - 1)*(rr - 1));
        rate = r2 + (r1 - r2)*qq;
    end
    
    % Reference m1 and m2
    m2ref = sigma00*rho(k);
    m1ref = rho(k) - m2ref;
    % Avoid negative transfers
    if m1(k) < 0.5*m1ref
        relabel.trialM21_tkeB(k) = 0;
        relabel.trialM12_tkeB(k) = m1(k)*rate;
%         dM21dm1(k) = 0;
%         dM21dm2(k) = 0;
%         dM12dm1(k) = rate;
%         dM12dm2(k) = 0;
    elseif m2(k) < 0.5*m2ref
        relabel.trialM12_tkeB(k) = 0;
        relabel.trialM21_tkeB(k) = m2(k)*rate;
%         dM21dm1(k) = 0;
%         dM21dm2(k) = rate;
%         dM12dm1(k) = 0;
%         dM12dm2(k) = 0;
    else
        relabel.trialM12_tkeB(k) = mf1(k)*(2*m2(k) - m2ref)*rate;
        relabel.trialM21_tkeB(k) = mf2(k)*(2*m1(k) - m1ref)*rate;
%         dM21dm1(k) = 2*sigma2(k)*rate;
%         dM21dm2(k) = (2*m1(k) - m1ref)*rate/rho(k);
%         dM12dm1(k) = (2*m2(k) - m2ref)*rate/rho(k);
%         dM12dm2(k) = 2*sigma1(k)*rate;
    end
    
end


% Interpolate entrainment and detrainment to w levels
% using `reversed' weighting for conservation
M12bar = weight_to_w(grid,M12);
M21bar = weight_to_w(grid,M21);

% Save computed values 
relabel.ideal = ones(1,nz)*sigma00;
relabel.M12 = M12;
relabel.M21 = M21;
relabel.M12bar = M12bar;
relabel.M21bar = M21bar;
relabel.dM21dm1 = dM21dm1;
relabel.dM21dm2 = dM21dm2;
relabel.dM21dw1 = dM21dw1;
relabel.dM21dw2 = dM21dw2;
relabel.dM21deta1 = dM21deta1;
relabel.dM21deta2 = dM21deta2;
relabel.dM21dq1 = dM21dq1;
relabel.dM21dq2 = dM21dq2;
relabel.dM12dm1 = dM12dm1;
relabel.dM12dm2 = dM12dm2;
relabel.dM12dw1 = dM12dw1;
relabel.dM12dw2 = dM12dw2;
relabel.dM12deta1 = dM12deta1;
relabel.dM12deta2 = dM12deta2;
relabel.dM12dq1 = dM12dq1;
relabel.dM12dq2 = dM12dq2;


% For testing
%disp('*** ZERO ENTRAINMENT ***')
%relabel.M21 = zeros(1,nz);
%relabel.M21bar = zeros(1,nzp);
%relabel.dM21dm1 = zeros(1,nz);
%relabel.dM21dm2 = zeros(1,nz);
%relabel.dM21dw1 = zeros(1,nz);
%relabel.dM21dw2 = zeros(1,nz);

%disp('*** ZERO DETRAINMENT ***')
%relabel.M12 = zeros(1,nz);
%relabel.M12bar = zeros(1,nzp);
%relabel.dM12dm1 = zeros(1,nz);
%relabel.dM12dm2 = zeros(1,nz);
%relabel.dM12dw1 = zeros(1,nz);
%relabel.dM12dw2 = zeros(1,nz);

% Entrained and detrained values of eta
[relabel.etahat12,relabel.detahat12deta1,relabel.detahat12deta2,...
 relabel.etahat21,relabel.detahat21deta1,relabel.detahat21deta2] ...
    = findqhat(eta1,eta2,bentraint,bdetraint);

% Entrained and detrained values of q
[relabel.qhat12,relabel.dqhat12dq1,relabel.dqhat12dq2,...
 relabel.qhat21,relabel.dqhat21dq1,relabel.dqhat21dq2] ...
    = findqhat(q1,q2,bentrainq,bdetrainq);

% Entrained and detrained values of w
[relabel.what12,relabel.dwhat12dw1,relabel.dwhat12dw2,...
 relabel.what21,relabel.dwhat21dw1,relabel.dwhat21dw2] ...
    = findqhat(w1,w2,bentrainw,bdetrainw);

% Entrained and detrained values of u and v
[relabel.uhat12,relabel.duhat12du1,relabel.duhat12du2,...
 relabel.uhat21,relabel.duhat21du1,relabel.duhat21du2] ...
    = findqhat(u1,u2,bentrainu,bdetrainu);
[relabel.vhat12,relabel.dvhat12dv1,relabel.dvhat12dv2,...
 relabel.vhat21,relabel.dvhat21dv1,relabel.dvhat21dv2] ...
    = findqhat(v1,v2,bentrainu,bdetrainu);


% ------------------

% Also compute some candidate entrainment and detrainment rates for
% comparison

% Preliminary calculations

% Inverse sqrt of 2
rroot2 = sqrt(0.5);

% sqrt of -1 * buoyancy frequency squared (where it's unstable)
% or buoyancy frequency (where its stable)
n1neg =     sqrt(max(-eos.nsq1,0));
n1pos = min(sqrt(max( eos.nsq1,0)),rdt);
n2neg =     sqrt(max(-eos.nsq2,0));
n2pos = min(sqrt(max( eos.nsq2,0)),rdt);

% Vertically averaged updraft w
w2bar = max(0,grid.abovep.*w2(2:nzp) + grid.belowp.*w2(1:nz));
% and vertical derivative
dw2dz = (w2(2:nzp) - w2(1:nz))./grid.dzp;

% Subfilter standard deviation of w, assuming tke has equal contributions
% from u, v and w
wstd = sqrt(tke2*2/3);
wstdw(2:nz) = grid.abovew(2:nz).*wstd(2:nz) + grid.beloww(2:nz).*wstd(1:nz-1);
wstdw(1) = wstdw(2);
wstdw(nzp) = wstdw(nz);
% And subfilter standard deviations of eta and q
etastd = sqrt(state.fluid(2).vareta);
qstd   = sqrt(state.fluid(2).varq);

% shear between fluids
du = u2 - u1;
dv = v2 - v1;
dww = w2 - w1;
dw = grid.abovep.*dww(2:nzp) + grid.belowp.*dww(1:nz);
absdu = sqrt(du.*du + dv.*dv + dw.*dw);

% Candidate relabelling rates

% Buoyancy-based formulations
relabel.trialM21_buoy = 0.5*m2.*max( buoybar,0)./(w2bar + 0.2*wstd);
relabel.trialM12_buoy = 0.5*m2.*max(-buoybar,0)./(w2bar + 0.2*wstd);

% Static stability-based formulations
relabel.trialM21_N1 = 0.2*m1.*n1neg;
relabel.trialM12_N1 = 0.2*m1.*n1pos;
relabel.trialM21_N2 = 0.2*m1.*n2neg;
relabel.trialM12_N2 = 0.2*m1.*n2pos;

% Mixing-based formulations
relabel.trialM21_du  = sigma1.*m2.*absdu./scales.L_turb2;
relabel.trialM12_du  = sigma1.*m2.*absdu./scales.L_turb2;
relabel.trialM21_tke = sigma1.*m2.*sqrt(tke2)./(scales.L_turb2);
relabel.trialM12_tke = sigma1.*m2.*sqrt(tke2)./(scales.L_turb2);
% relabel.trialM21_tkeB = sigma1.*m2.*sqrt(tke2)./(scales.L_plume);
% relabel.trialM12_tkeB = sigma1.*m2.*sqrt(tke2)./(scales.L_plume);
relabel.trialM21_K = sigma1.*m2.*kdiffw2./(scales.L_turb2.^2);
relabel.trialM12_K = sigma1.*m2.*kdiffw2./(scales.L_turb2.^2);

% APDF-based formulation
relabel.w2bar = w2bar;
relabel.wstd = wstd;
% 1/sqrt(2) * normalized variable at cut of pdf on p-levels
chi_cut = - rroot2*w2bar./wstd;
% and w-levels (bounded to prevent divide by zero below)
chi_cut_w = max(-rroot2*w2./wstdw,-4);

% Fraction of pdf cut off
%relabel.frac = 0.5.*(1 + erf(-rroot2*w2bar./wstd));
relabel.frac = 0.5*(1 + erf(chi_cut));

% Rate at which it is cut off
%rate_apdf = 5*(wstd./scales.L_plume);
rate_apdf = min(10*max(0,-dw2dz),rdt);

% Resulting detrainment rate
relabel.trialM12_apdf = rate_apdf.*m2.*relabel.frac;

% Trial formulation - a combination of the above
relabel.trialM21 = relabel.trialM21_N1   + relabel.trialM21_tkeB;
relabel.trialM12 = relabel.trialM12_apdf + relabel.trialM12_tkeB;

% w, eta and q of detrained fluid
% wbyw = min(rroot2*w2./wstdw,4);  % Bound by 4 to avoid divide by zero in next line
% relabel.trialwhat12 = w2 - sqrt(2/pi)*wstdw.*exp(-wbyw.^2)./(1 + erf(-wbyw));
chi_hat = - sqrt(2/pi)*exp(-chi_cut_w.^2)./(1 + erf(chi_cut_w));
relabel.trialwhat12_apdf   = w2   + wstdw .*chi_hat;
relabel.trialetahat12_apdf = eta2 + etastd.*chi_hat;
relabel.trialqhat12_apdf   = q2   + qstd  .*chi_hat;


% Try candidate entrainment/detrainment
relabel.M21 = relabel.trialM21;
relabel.M12 = relabel.trialM12;
% relabel.what12   = relabel.trialwhat12_apdf;
% relabel.etahat12 = relabel.trialetahat12_apdf;
% relabel.qhat12   = relabel.trialqhat12_apdf;


% relabel.M21 = zeros(1,nz);
% relabel.dM21dm1 = zeros(1,nz);
% relabel.dM21dm2 = zeros(1,nz);
% relabel.dM21dw1 = zeros(1,nz);
% relabel.dM21dw2 = zeros(1,nz);
% relabel.dM21deta1 = zeros(1,nz);
% relabel.dM21deta2 = zeros(1,nz);
% relabel.dM21dq1 = zeros(1,nz);
% relabel.dM21dq2 = zeros(1,nz);

% Try candidate detrainment
%relabel.M12 = relabel.trialM12s + relabel.trialM12m;
%relabel.dM12dm1 = zeros(1,nz);
%relabel.dM12dm2 = relabel.M12./m2;
%relabel.dM12dw1 = zeros(1,nz);
%relabel.dM12dw2 = zeros(1,nz);
%relabel.dM12deta1 = zeros(1,nz);
%relabel.dM12deta2 = zeros(1,nz);
%relabel.dM12dq1 = zeros(1,nz);
%relabel.dM12dq2 = zeros(1,nz);

% -----------------------------------------------------

% Interpolate entrainment and detrainment to w levels
% using `reversed' weighting for conservation
relabel.M12bar = weight_to_w(grid,relabel.M12);
relabel.M21bar = weight_to_w(grid,relabel.M21);



end