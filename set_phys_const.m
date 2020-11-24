function phys = set_phys_const( )

% Set physical constants

phys.gravity  = 9.806;     % Gravity
phys.k        = 0.4;       % von Karman constant
phys.zrough   = 0.035;     % Roughness length for diffusion
phys.p00      = 100000;    % Reference pressure for potential temperature
phys.coriolis = 8.5e-5;    % Coriolis parameter

end

