% Estimate profiles of fluid 1 and 2 cloud fraction,
% Hence estimate cloud base and top and cloud cover


% Stand-alone method - used for checking or if find_eos is used rather
% than find_eos_sg.
%
% % Assume etastd and qstd have been computed already
% 
% % Assumed correlation between entropy and total water
% % *** Set this in constants.params ***
% rr = 0.0;
% 
% % Useful constant
% rr2 = sqrt(0.5);
% 
% % Loop over levels
% for k = 1:nzp
%     
%     % w level pressure
%     if k == 1
%         pbar = grid.extrapb1*p(1) + grid.extrapb2*p(2);
%     elseif k == nzp
%         pbar = grid.extraptnz*p(nz) + grid.extraptnzm*p(nz-1);
%     else
%         pbar   = grid.abovew(k)*p(k) ...
%                + grid.beloww(k)*p(k-1);
%     end
%     
%     % Fluid 1
%     % Compute qsat and dqsat/deta at this eta
%     [qsat,dqsatdeta] = find_qsatl(pbar,eta1(k),Tw1(k),constants.therm);
%     
%     % Difference between mean q and qsat
%     Deltaq = q1(k) - qsat;
%     
%     % Standard deviation parameter
%     etastd = sqrt(state_new.fluid(1).vareta(k));
%     qstd   = sqrt(state_new.fluid(1).varq(k));
%     sq2 = qstd*qstd ...
%         - 2*rr*qstd*etastd*dqsatdeta ...
%         + etastd*etastd*dqsatdeta*dqsatdeta;
%     sq = sqrt(sq2);
%     
%     % Q1 parameter
%     Q1 = Deltaq/sq;
%     
%     % Cloud fraction
%     cld_frac1x(k) = sigma1w(k)*0.5*(1 + erf(rr2*Q1));
%     
%     % Fluid 2
%     % Compute qsat and dqsat/deta at this eta
%     [qsat,dqsatdeta] = find_qsatl(pbar,eta2(k),Tw2(k),constants.therm);
%     
%     % Difference between mean q and qsat
%     Deltaq = q2(k) - qsat;
%     
%     % Standard deviation parameter
%     etastd = sqrt(state_new.fluid(2).vareta(k));
%     qstd   = sqrt(state_new.fluid(2).varq(k));
%     sq2 = qstd*qstd ...
%         - 2*rr*qstd*etastd*dqsatdeta ...
%         + etastd*etastd*dqsatdeta*dqsatdeta;
%     sq = sqrt(sq2);
%     
%     % Q1 parameter
%     Q1 = Deltaq/sq;
%     
%     % Cloud fraction
%     cld_frac2x(k) = sigma2w(k)*0.5*(1 + erf(rr2*Q1));
%     
% end


% If we're using find_eos_sg cloud fractions are already diagnosed
cld_frac1 = sigma1w.*eos.cldfrac1;
cld_frac2 = sigma2w.*eos.cldfrac2;


tot_cld_frac = cld_frac1 + cld_frac2;

% Total cloud cover
tot_cld_cov = max(tot_cld_frac);

% Threshold for defining base and top
% cld_thresh = 0.01;
cld_thresh = 0.001;
z_cld_top = 0;
z_cld_base = 0;
for k = 1:nzp
    if tot_cld_frac(k) > cld_thresh
        z_cld_top = grid.zw(k);
    end
    rk = nz+2-k;
    if tot_cld_frac(rk) > cld_thresh
        z_cld_base = grid.zw(rk);
    end
end
% for k = 1:nzp
    % if eos.cldfrac1(k) > cld_thresh | eos.cldfrac2(k) > cld_thresh
        % z_cld_top = grid.zw(k);
    % end
    % rk = nz+2-k;
    % if eos.cldfrac1(rk) > cld_thresh | eos.cldfrac2(rk) > cld_thresh
        % z_cld_base = grid.zw(rk);
    % end
% end

