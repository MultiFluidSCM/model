% Integrate the governing equations in time

for istep = 1:time.nstop
    
    time.istep = istep;
    
    % Step forward
    advance
    time.t = time.t + time.dt;
    
    % Compute diagnostics and produce plots
    if mod(istep,20) == 0
        gdiags = global_diags(grid,state_new,constants);
        plottype = 0;
        plot_diagnostics
        figure(1)
    end
    
    % Save a restart file if needed
    if mod(istep,60000) == 0
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