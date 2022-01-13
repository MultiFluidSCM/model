% Update timeseries for cloud properties
if exist('ts')
    diagnose_cloud_frac
    
    index = 1;
    if isfield(ts, 'time_high_res')
        index = length(ts.time_high_res) + 1;
    end
    
    ts.time_high_res(index) = time.t;
    ts.zstar(index) = scales.zstar;
    ts.zcbaseSG(index) = z_cld_base;
    ts.zctopSG(index) = z_cld_top;
    ts.totcldcov(index) = tot_cld_cov;
    ts.cloud_fraction(:,index) = tot_cld_frac;
    ts.cloud_fraction1(:,index) = eos.cldfrac1;
    ts.cloud_fraction2(:,index) = eos.cldfrac2;
    ts.cloud_fraction1_sigma1(:,index) = eos.cldfrac1.*sigma1w;
    ts.cloud_fraction2_sigma2(:,index) = eos.cldfrac2.*sigma2w;
    
end
