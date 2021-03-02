% Update timeseries for cloud properties
if exist('ts')
    diagnose_cloud_frac
    
    index = round(time.t/time.dt);
    ts.time_high_res(index) = time.t;
    ts.zstar(index) = scales.zstar;
    ts.zcbaseSG(index) = z_cld_base;
    ts.zctopSG(index) = z_cld_top;
    ts.totcldcov(index) = tot_cld_cov;
end