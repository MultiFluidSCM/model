% Integrate the governing equations in time

for istep = 1:time.nstop
    
    time.istep = istep;
    
    % Step forward if simulation hasn't crashed
    if ~isnan(state_new.fluid(1).m(1))
        advance
    end
    time.t = time.t + time.dt;
    
    gdiags = global_diags(grid,state_new,constants);
    SCM_energy_time(istep)    = time.t;
    SCM_energy_total(istep)   = gdiags.energy1 + gdiags.energy2;
    SCM_energy_surface(istep) = time.dt*(force.sshf + 2.5e6*force.sqf);
    
    % Compute diagnostics and produce plots
    if mod(istep,20) == 0
        % gdiags = global_diags(grid,state_new,constants);
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