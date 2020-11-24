function check = check_eos(grid, state, constants)

% Determine densities and find residuals in equations of state

% Also compute terms needed in linearization. The quantity needed
% is the derivative of (rho' / rho^2) wrt p', eta', and q'

nz = grid.nz;
nzp = nz + 1;

% On p levels
for k = 1
    
    % Fluid 1
    qbar   = grid.abovep(k)*state.fluid(1).q(k+1) ...
           + grid.belowp(k)*state.fluid(1).q(k);
    etabar = grid.abovep(k)*state.fluid(1).eta(k+1) ...
           + grid.belowp(k)*state.fluid(1).eta(k);
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(state.p(k),          ...
                                                   state.fluid(1).T(k), ...
                                                   qbar,                ...
                                                   constants.therm);
    check.p.p = state.p(k);    
    check.p1.eta = etabar;
    check.p1.T = state.fluid(1).T(k);
    check.p1.q = qbar;
    check.p1.g = g;
    check.p1.gw = gw;
    check.p1.alpha = gp;
    check.p1.gtt = gtt;
    check.p1.gtw = gtw;
    check.p1.gpt = gpt;
    check.p1.drdeta = gpt/gtt;
    check.p1.drdq = (gpt*gtw/gtt - gpw);
    check.p1.drdp = (gpt*gpt/gtt - gpp);
    check.p1.gpw = gpw;
    check.p1.gpp = gpp;
    
    % Fluid 2
    qbar   = grid.abovep(k)*state.fluid(2).q(k+1) ...
           + grid.belowp(k)*state.fluid(2).q(k);
    etabar = grid.abovep(k)*state.fluid(2).eta(k+1) ...
           + grid.belowp(k)*state.fluid(2).eta(k);
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(state.p(k),          ...
                                                   state.fluid(2).T(k), ...
                                                   qbar,                ...
                                                   constants.therm);
    check.p2.eta = etabar;
    check.p2.T = state.fluid(2).T(k);
    check.p2.q = qbar;
    check.p2.g = g;
    check.p2.gw = gw;
    check.p2.alpha = gp;
    check.p2.gtt = gtt;
    check.p2.gtw = gtw;
    check.p2.gpt = gpt;
    check.p2.drdeta = gpt/gtt;
    check.p2.drdq = (gpt*gtw/gtt - gpw);
    check.p2.drdp = (gpt*gpt/gtt - gpp);
    check.p2.gpw = gpw;
    check.p2.gpp = gpp;

end


% On w levels
for k = 1
    
    if k == 1
        pbar = grid.extrapb1*state.p(1) + grid.extrapb2*state.p(2);
    elseif k == nzp
        pbar = grid.extraptnz*state.p(nz) + grid.extraptnzm*state.p(nz-1);
    else
        pbar   = grid.abovew(k)*state.p(k) ...
               + grid.beloww(k)*state.p(k-1);
    end
    
    % Fluid 1
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(pbar,                 ...
                                                   state.fluid(1).Tw(k), ...
                                                   state.fluid(1).q(k),  ...
                                                   constants.therm);
    check.w.p = pbar;    
    check.w1.eta = state.fluid(1).eta(k);
    check.w1.T = state.fluid(1).Tw(k);
    check.w1.q = state.fluid(1).q(k);
    check.w1.g = g;
    check.w1.gw = gw;
    check.w1.alpha = gp;
    check.w1.gtt = gtt;
    check.w1.gtw = gtw;
    check.w1.gpt = gpt;
    check.w1.drdeta = gpt/gtt;
    check.w1.drdq = (gpt*gtw/gtt - gpw);
    check.w1.drdp = (gpt*gpt/gtt - gpp);
    check.w1.gpw = gpw;
    check.w1.gpp = gpp;
    
    % Fluid2
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(pbar,                 ...
                                                   state.fluid(2).Tw(k), ...
                                                   state.fluid(2).q(k),  ...
                                                   constants.therm);
    check.w2.eta = state.fluid(2).eta(k);
    check.w2.T = state.fluid(2).Tw(k);
    check.w2.q = state.fluid(2).q(k);
    check.w2.g = g;
    check.w2.gw = gw;
    check.w2.alpha = gp;
    check.w2.gtt = gtt;
    check.w2.gtw = gtw;
    check.w2.gpt = gpt;
    check.w2.drdeta = gpt/gtt;
    check.w2.drdq = (gpt*gtw/gtt - gpw);
    check.w2.drdp = (gpt*gpt/gtt - gpp);
    check.w2.gpw = gpw;
    check.w2.gpp = gpp;
    
end
