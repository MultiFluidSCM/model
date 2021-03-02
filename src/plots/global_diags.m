function [ gdiags ] = global_diags( grid, state, constants )
%GLOBAL_DIAGS Compute a variety of global integral
% diagnostics

m1bar = weight_to_w(grid,state.fluid(1).m);
m2bar = weight_to_w(grid,state.fluid(2).m);

% Mass
mass1 = sum(state.fluid(1).m.*grid.dzp);
mass2 = sum(state.fluid(2).m.*grid.dzp);

% Entropy
entropy1 = sum(m1bar.*state.fluid(1).eta.*grid.dzw);
entropy2 = sum(m2bar.*state.fluid(2).eta.*grid.dzw);

% Water
water1 = sum(m1bar.*state.fluid(1).q.*grid.dzw);
water2 = sum(m2bar.*state.fluid(2).q.*grid.dzw);

% Energy
internal1  = 0;
potential1 = 0;
kinetic1   = 0;
turbKE1    = 0;
internal2  = 0;
potential2 = 0;
kinetic2   = 0;
turbKE2    = 0;
for k = 1:grid.nz
    pp = state.p(k);
    tt = state.fluid(1).T(k);
    dm = state.fluid(1).m(k)*grid.dzp(k);
    qbar = grid.aboves(k)*state.fluid(1).q(k+1) + grid.belows(k)*state.fluid(1).q(k);
    wsqbar = grid.aboves(k)*state.fluid(1).w(k+1)^2 ...
           + grid.belows(k)*state.fluid(1).w(k  )^2;
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = ...
            gibbs(pp,tt,qbar,constants.therm);
    internal1  = internal1  + (g - pp*gp - tt*gt)*dm;
    potential1 = potential1 + constants.phi(k)*dm;
    kinetic1   = kinetic1   + 0.5*(wsqbar + state.fluid(1).u(k)^2 + state.fluid(1).v(k)^2)*dm;
    turbKE1    = turbKE1    + state.fluid(1).tke(k)*dm;
    tt = state.fluid(2).T(k);
    dm = state.fluid(2).m(k)*grid.dzp(k);
    qbar = grid.aboves(k)*state.fluid(2).q(k+1) + grid.belows(k)*state.fluid(2).q(k);
    wsqbar = grid.aboves(k)*state.fluid(2).w(k+1)^2 ...
           + grid.belows(k)*state.fluid(2).w(k  )^2;
    [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = ...
            gibbs(pp,tt,qbar,constants.therm);
    internal2  = internal2  + (g - pp*gp - tt*gt)*dm;
    potential2 = potential2 + constants.phi(k)*dm;
    kinetic2   = kinetic2   + 0.5*(wsqbar + state.fluid(2).u(k)^2 + state.fluid(2).v(k)^2)*dm;
    turbKE2    = turbKE2    + state.fluid(2).tke(k)*dm;
end
energy1 = internal1 + potential1 + kinetic1 + turbKE1;
energy2 = internal2 + potential2 + kinetic2 + turbKE2;

% Display diagnostics
disp([' '])
disp(['Diagnostics:'])
disp(['Variable      Fluid 1   Fluid 2   Total'])
disp(['Mass         ' num2str(mass1     ) '  ' num2str(mass2     ) '  ' num2str(mass1     + mass2     )])
disp(['Entropy      ' num2str(entropy1  ) '  ' num2str(entropy2  ) '  ' num2str(entropy1  + entropy2  )])
disp(['Water        ' num2str(water1    ) '  ' num2str(water2    ) '  ' num2str(water1    + water2    )])
disp(['Energy(I)    ' num2str(internal1 ) '  ' num2str(internal2 ) '  ' num2str(internal1 + internal2 )])
disp(['Energy(P)    ' num2str(potential1) '  ' num2str(potential2) '  ' num2str(potential1+ potential2)])
disp(['Energy(K)    ' num2str(kinetic1  ) '  ' num2str(kinetic2  ) '  ' num2str(kinetic1  + kinetic2  )])
disp(['Energy(TKE)  ' num2str(turbKE1   ) '  ' num2str(turbKE2   ) '  ' num2str(turbKE1   + turbKE2   )])
disp(['Energy(Tot)  ' num2str(energy1   ) '  ' num2str(energy2   ) '  ' num2str(energy1   + energy2   )])
disp([' '])

% Save diagnostics for use elsewhere
gdiags.mass1 = mass1;
gdiags.mass2 = mass2;
gdiags.entropy1 = entropy1;
gdiags.entropy2 = entropy2;
gdiags.water1 = water1;
gdiags.water2 = water2;
gdiags.internal1 = internal1;
gdiags.internal2 = internal2;
gdiags.potential1 = potential1;
gdiags.potential2 = potential2;
gdiags.kinetic1 = kinetic1;
gdiags.kinetic2 = kinetic2;
gdiags.tke1 = turbKE1;
gdiags.tke2 = turbKE2;
gdiags.energy1 = energy1;
gdiags.energy2 = energy2;

end
