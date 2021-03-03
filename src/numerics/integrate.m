% Integrate the governing equations in time

for istep = 1:time.nstop
    
    time.istep = istep;
    
    % Step forward
    advance
    time.t = time.t + time.dt;
    
    % Compute diagnostics and produce plots
    update_timeseries_cloud
    if mod(istep,60) == 0
        gdiags = global_diags(grid,state_new,constants);
        plottype = 0;
        plot_diagnostics
    end
    
    % Save a restart file if needed
    if mod(istep,10800) == 0
        save_restart
    end
    
    % Save data for comparison with LES
    compare_LES
    
    % Pause to look at diagnostics
    if mod(istep,1000000) == 0
        pause
    end
    
end

% Final plots
plottype = 0;
plot_diagnostics