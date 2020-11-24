function constants = set_constants(grid)

% Set physical, thermodynamic, and paramaterization constants

constants.phys  = set_phys_const();
constants.therm = set_therm_const();
constants.param = set_param_const();

% Compute profile of geopotential
constants.phi = constants.phys.gravity * grid.zp;

end

