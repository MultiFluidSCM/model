% Estimate profiles of fluid 1 and 2 cloud fraction,
% Hence estimate cloud base and top and cloud cover


% If we're using find_eos_sg cloud fractions are already diagnosed
cld_frac1 = sigma1w.*eos.cldfrac1;
cld_frac2 = sigma2w.*eos.cldfrac2;


tot_cld_frac = cld_frac1 + cld_frac2;

% Total cloud cover
tot_cld_cov = max(tot_cld_frac);

% Threshold for defining base and top
cld_thresh = settings.constants.param.cld_thresh;
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

