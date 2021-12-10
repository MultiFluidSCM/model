% Check that the implied tendencies of specific quantities are
% equal when the properties of the two fluids are equal, even
% if sigm1 and sigma2 depend on z.

% The discretization should ensure all contributions are equal
% except for the transport of w-level tracers; this latter
% cannot be achieved while retaining a consistent w-level mass
% budget, but the effects should be minor.

% Call equate_fluid_properties before calling tendencies,
% then call this routine from the end of tendencies

kkk = 1:4;

% Entropy
disp(' ')
disp('diff etadot')
disp('transport')
m1bardot = weight_to_w(grid,tend.fluid(1).m.transport);
m2bardot = weight_to_w(grid,tend.fluid(2).m.transport);
xxx1 = (tend.fluid(1).meta.transport - eta1.*m1bardot)./m1bar;
xxx2 = (tend.fluid(2).meta.transport - eta2.*m2bardot)./m2bar;
xxx2(kkk) - xxx1(kkk)
disp('relabel')
m1bardot = weight_to_w(grid,tend.fluid(1).m.relabel);
m2bardot = weight_to_w(grid,tend.fluid(2).m.relabel);
xxx1 = (tend.fluid(1).meta.relabel - eta1.*m1bardot)./m1bar;
xxx2 = (tend.fluid(2).meta.relabel - eta2.*m2bardot)./m2bar;
xxx2(kkk) - xxx1(kkk)
disp('diffuse')
xxx1 = tend.fluid(1).meta.diffuse./m1bar;
xxx2 = tend.fluid(2).meta.diffuse./m2bar;
xxx2(kkk) - xxx1(kkk)
disp('diffent')
xxx1 = tend.fluid(1).meta.diffent./m1bar;
xxx2 = tend.fluid(2).meta.diffent./m2bar;
xxx2(kkk) - xxx1(kkk)
disp('dissn')
xxx1 = tend.fluid(1).meta.dissn./m1bar;
xxx2 = tend.fluid(2).meta.dissn./m2bar;
xxx2(kkk) - xxx1(kkk)
disp('tot')
m1bardot = weight_to_w(grid,tend.fluid(1).m.tot);
m2bardot = weight_to_w(grid,tend.fluid(2).m.tot);
xxx1 = (tend.fluid(1).meta.tot - eta1.*m1bardot)./m1bar;
xxx2 = (tend.fluid(2).meta.tot - eta2.*m2bardot)./m2bar;
xxx2(kkk) - xxx1(kkk)

% Water
disp(' ')
disp('diff qdot')
disp('transport')
m1bardot = weight_to_w(grid,tend.fluid(1).m.transport);
m2bardot = weight_to_w(grid,tend.fluid(2).m.transport);
xxx1 = (tend.fluid(1).mq.transport - q1.*m1bardot)./m1bar;
xxx2 = (tend.fluid(2).mq.transport - q2.*m2bardot)./m2bar;
xxx2(kkk) - xxx1(kkk)
disp('relabel')
m1bardot = weight_to_w(grid,tend.fluid(1).m.relabel);
m2bardot = weight_to_w(grid,tend.fluid(2).m.relabel);
xxx1 = (tend.fluid(1).mq.relabel - q1.*m1bardot)./m1bar;
xxx2 = (tend.fluid(2).mq.relabel - q2.*m2bardot)./m2bar;
xxx2(kkk) - xxx1(kkk)
disp('diffuse')
xxx1 = tend.fluid(1).mq.diffuse./m1bar;
xxx2 = tend.fluid(2).mq.diffuse./m2bar;
xxx2(kkk) - xxx1(kkk)
disp('diffent')
xxx1 = tend.fluid(1).mq.diffent./m1bar;
xxx2 = tend.fluid(2).mq.diffent./m2bar;
xxx2(kkk) - xxx1(kkk)
disp('tot')
m1bardot = weight_to_w(grid,tend.fluid(1).m.tot);
m2bardot = weight_to_w(grid,tend.fluid(2).m.tot);
xxx1 = (tend.fluid(1).mq.tot - q1.*m1bardot)./m1bar;
xxx2 = (tend.fluid(2).mq.tot - q2.*m2bardot)./m2bar;
xxx2(kkk) - xxx1(kkk)

% Vertical velocity
disp(' ')
disp('diff wdot')
disp('transport')
m1bardot = weight_to_w(grid,tend.fluid(1).m.transport);
m2bardot = weight_to_w(grid,tend.fluid(2).m.transport);
xxx1 = (tend.fluid(1).mw.transport - w1.*m1bardot)./m1bar;
xxx2 = (tend.fluid(2).mw.transport - w2.*m2bardot)./m2bar;
xxx2(kkk) - xxx1(kkk)
disp('relabel')
m1bardot = weight_to_w(grid,tend.fluid(1).m.relabel);
m2bardot = weight_to_w(grid,tend.fluid(2).m.relabel);
xxx1 = (tend.fluid(1).mw.relabel - w1.*m1bardot)./m1bar;
xxx2 = (tend.fluid(2).mw.relabel - w2.*m2bardot)./m2bar;
xxx2(kkk) - xxx1(kkk)
disp('diffuse')
xxx1 = tend.fluid(1).mw.diffuse./m1bar;
xxx2 = tend.fluid(2).mw.diffuse./m2bar;
xxx2(kkk) - xxx1(kkk)
disp('diffent')
xxx1 = tend.fluid(1).mw.diffent./m1bar;
xxx2 = tend.fluid(2).mw.diffent./m2bar;
xxx2(kkk) - xxx1(kkk)
disp('pgterm')
xxx1 = tend.fluid(1).mw.pgterm./m1bar;
xxx2 = tend.fluid(2).mw.pgterm./m2bar;
% xxx2(kkk)
% xxx1(kkk)
xxx2(kkk) - xxx1(kkk)
disp('drag')
xxx1 = tend.fluid(1).mw.drag./m1bar;
xxx2 = tend.fluid(2).mw.drag./m2bar;
xxx2(kkk) - xxx1(kkk)
disp('tot')
m1bardot = weight_to_w(grid,tend.fluid(1).m.tot);
m2bardot = weight_to_w(grid,tend.fluid(2).m.tot);
xxx1 = (tend.fluid(1).mw.tot - w1.*m1bardot)./m1bar;
xxx2 = (tend.fluid(2).mw.tot - w2.*m2bardot)./m2bar;
xxx2(kkk) - xxx1(kkk)

% Horizontal velocity
disp(' ')
disp('diff udot')
disp('transport')
m1dot = tend.fluid(1).m.transport;
m2dot = tend.fluid(2).m.transport;
xxx1 = (tend.fluid(1).mu.transport - u1.*m1dot)./m1;
xxx2 = (tend.fluid(2).mu.transport - u2.*m2dot)./m2;
xxx2(kkk) - xxx1(kkk)
disp('relabel')
m1dot = tend.fluid(1).m.relabel;
m2dot = tend.fluid(2).m.relabel;
xxx1 = (tend.fluid(1).mu.relabel - u1.*m1dot)./m1;
xxx2 = (tend.fluid(2).mu.relabel - u2.*m2dot)./m2;
xxx2(kkk) - xxx1(kkk)
disp('diffuse')
xxx1 = tend.fluid(1).mu.diffuse./m1;
xxx2 = tend.fluid(2).mu.diffuse./m2;
xxx2(kkk) - xxx1(kkk)
disp('diffent')
xxx1 = tend.fluid(1).mu.diffent./m1;
xxx2 = tend.fluid(2).mu.diffent./m2;
xxx2(kkk) - xxx1(kkk)
disp('coriolis')
xxx1 = tend.fluid(1).mu.coriolis./m1;
xxx2 = tend.fluid(2).mu.coriolis./m2;
xxx2(kkk) - xxx1(kkk)
disp('tot')
m1bardot = tend.fluid(1).m.tot;
m2bardot = tend.fluid(2).m.tot;
xxx1 = (tend.fluid(1).mu.tot - u1.*m1dot)./m1;
xxx2 = (tend.fluid(2).mu.tot - u2.*m2dot)./m2;
xxx2(kkk) - xxx1(kkk)

disp(' ')
disp('diff vdot')
disp('transport')
m1dot = tend.fluid(1).m.transport;
m2dot = tend.fluid(2).m.transport;
xxx1 = (tend.fluid(1).mv.transport - v1.*m1dot)./m1;
xxx2 = (tend.fluid(2).mv.transport - v2.*m2dot)./m2;
xxx2(kkk) - xxx1(kkk)
disp('relabel')
m1dot = tend.fluid(1).m.relabel;
m2dot = tend.fluid(2).m.relabel;
xxx1 = (tend.fluid(1).mv.relabel - v1.*m1dot)./m1;
xxx2 = (tend.fluid(2).mv.relabel - v2.*m2dot)./m2;
xxx2(kkk) - xxx1(kkk)
disp('diffuse')
xxx1 = tend.fluid(1).mv.diffuse./m1;
xxx2 = tend.fluid(2).mv.diffuse./m2;
xxx2(kkk) - xxx1(kkk)
disp('diffent')
xxx1 = tend.fluid(1).mv.diffent./m1;
xxx2 = tend.fluid(2).mv.diffent./m2;
xxx2(kkk) - xxx1(kkk)
disp('coriolis')
xxx1 = tend.fluid(1).mv.coriolis./m1;
xxx2 = tend.fluid(2).mv.coriolis./m2;
xxx2(kkk) - xxx1(kkk)
disp('tot')
m1bardot = tend.fluid(1).m.tot;
m2bardot = tend.fluid(2).m.tot;
xxx1 = (tend.fluid(1).mv.tot - v1.*m1dot)./m1;
xxx2 = (tend.fluid(2).mv.tot - v2.*m2dot)./m2;
xxx2(kkk) - xxx1(kkk)

% TKE
disp(' ')
disp('diff tkedot')
disp('transport')
m1dot = tend.fluid(1).m.transport;
m2dot = tend.fluid(2).m.transport;
xxx1 = (tend.fluid(1).mtke.transport - tke1.*m1dot)./m1;
xxx2 = (tend.fluid(2).mtke.transport - tke2.*m2dot)./m2;
xxx2(kkk) - xxx1(kkk)
disp('relabel')
m1dot = tend.fluid(1).m.relabel;
m2dot = tend.fluid(2).m.relabel;
xxx1 = (tend.fluid(1).mtke.relabel - tke1.*m1dot)./m1;
xxx2 = (tend.fluid(2).mtke.relabel - tke2.*m2dot)./m2;
xxx2(kkk) - xxx1(kkk)
disp('diffuse')
xxx1 = tend.fluid(1).mtke.diffuse./m1;
xxx2 = tend.fluid(2).mtke.diffuse./m2;
xxx2(kkk) - xxx1(kkk)
disp('diffent')
xxx1 = tend.fluid(1).mtke.diffent./m1;
xxx2 = tend.fluid(2).mtke.diffent./m2;
xxx2(kkk) - xxx1(kkk)
disp('shear')
xxx1 = tend.fluid(1).mtke.shear./m1;
xxx2 = tend.fluid(2).mtke.shear./m2;
xxx2(kkk) - xxx1(kkk)
disp('bflux')
xxx1 = tend.fluid(1).mtke.bflux./m1;
xxx2 = tend.fluid(2).mtke.bflux./m2;
xxx2(kkk) - xxx1(kkk)
disp('drag')
xxx1 = tend.fluid(1).mtke.drag./m1;
xxx2 = tend.fluid(2).mtke.drag./m2;
xxx2(kkk) - xxx1(kkk)
disp('tot')
m1bardot = tend.fluid(1).m.tot;
m2bardot = tend.fluid(2).m.tot;
xxx1 = (tend.fluid(1).mtke.tot - tke1.*m1dot)./m1;
xxx2 = (tend.fluid(2).mtke.tot - tke2.*m2dot)./m2;
xxx2(kkk) - xxx1(kkk)


% eta variance (advective form)
disp(' ')
disp('diff varetadot')
disp('relabel')
xxx1 = tend.fluid(1).mvareta.relabel./m1;
xxx2 = tend.fluid(2).mvareta.relabel./m2;
xxx2(kkk) - xxx1(kkk)
disp('diffuse')
xxx1 = tend.fluid(1).mvareta.diffuse./m1;
xxx2 = tend.fluid(2).mvareta.diffuse./m2;
xxx2(kkk) - xxx1(kkk)
disp('diffent')
xxx1 = tend.fluid(1).mvareta.diffent./m1;
xxx2 = tend.fluid(2).mvareta.diffent./m2;
xxx2(kkk) - xxx1(kkk)
disp('dissn')
xxx1 = tend.fluid(1).mvareta.dissn./m1;
xxx2 = tend.fluid(2).mvareta.dissn./m2;
xxx2(kkk) - xxx1(kkk)
disp('tot')
xxx1 = tend.fluid(1).mvareta.tot./m1;
xxx2 = tend.fluid(2).mvareta.tot./m2;
xxx2(kkk) - xxx1(kkk)

% q variance (advective form)
disp(' ')
disp('diff varqdot')
disp('relabel')
xxx1 = tend.fluid(1).mvarq.relabel./m1;
xxx2 = tend.fluid(2).mvarq.relabel./m2;
xxx2(kkk) - xxx1(kkk)
disp('diffuse')
xxx1 = tend.fluid(1).mvarq.diffuse./m1;
xxx2 = tend.fluid(2).mvarq.diffuse./m2;
xxx2(kkk) - xxx1(kkk)
disp('diffent')
xxx1 = tend.fluid(1).mvarq.diffent./m1;
xxx2 = tend.fluid(2).mvarq.diffent./m2;
xxx2(kkk) - xxx1(kkk)
disp('dissn')
xxx1 = tend.fluid(1).mvarq.dissn./m1;
xxx2 = tend.fluid(2).mvarq.dissn./m2;
xxx2(kkk) - xxx1(kkk)
disp('tot')
xxx1 = tend.fluid(1).mvarq.tot./m1;
xxx2 = tend.fluid(2).mvarq.tot./m2;
xxx2(kkk) - xxx1(kkk)
