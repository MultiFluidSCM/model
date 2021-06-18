% Update timeseries for cloud properties
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