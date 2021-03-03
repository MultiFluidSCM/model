% Disect residuals into contributions from different tendencies

rr_transport = adt*(tend.fluid(2).mw.transport - tend_tru.fluid(2).mw.transport);
rr_diffuse   = adt*(tend.fluid(2).mw.diffuse   - tend_tru.fluid(2).mw.diffuse  );
rr_drag      = adt*(tend.fluid(2).mw.drag      - tend_tru.fluid(2).mw.drag     );
rr_pgterm    = adt*(tend.fluid(2).mw.pgterm    - tend_tru.fluid(2).mw.pgterm   );
rr_relabel   = adt*(tend.fluid(2).mw.relabel   - tend_tru.fluid(2).mw.relabel  );
x2bar = weight_to_w(grid,res2m);
rr_mbar      = -work.w2ubar.*x2bar;
rr_td = adt*m2bar.*dpdz.*eos.drdeta2.*eos.res_eta2;

rr = rr_transport + rr_diffuse + rr_drag + rr_pgterm + rr_relabel + rr_mbar + rr_td;

figure(12)
subplot(2,3,1)
plot(rr_transport(prange)+rr_mbar(prange),grid.zw(prange),'r')
title('transport + mbar')
set(gca,'fontsize',fs)
subplot(2,3,2)
plot(rr_drag(prange),grid.zw(prange),'r',rr_diffuse(prange),grid.zw(prange),'r--')
title('drag + diffuse')
set(gca,'fontsize',fs)
%subplot(2,3,3)
%plot(rr_td(prange),grid.zw(prange),'r',rr_mbar(prange),grid.zw(prange),'r--')
%title('thermo + mbar')
%set(gca,'fontsize',fs)
subplot(2,3,4)
plot(rr_pgterm(prange)+rr_td(prange),grid.zw(prange),'r',...
     rr_pgterm(prange),grid.zw(prange),'b',...
     rr_td(prange),grid.zw(prange),'g')
title('pgterm')
set(gca,'fontsize',fs)
subplot(2,3,3)
plot(rr_relabel(prange),grid.zw(prange),'r')
title('relabel')
set(gca,'fontsize',fs)
subplot(2,3,5)
plot(rr_transport(prange)+rr_mbar(prange)+rr_relabel(prange),grid.zw(prange),'r')
title('trans+mbar+relab')
set(gca,'fontsize',fs)
subplot(2,3,6)
plot(rr(prange),grid.zw(prange),'r')
title('total')
set(gca,'fontsize',fs)
%pause
%figure(1)




rr_transport = adt*(tend.fluid(2).m.transport - tend_tru.fluid(2).m.transport);
rr_relabel   = adt*(tend.fluid(2).m.relabel   - tend_tru.fluid(2).m.relabel  );

figure(13)
subplot(1,3,1)
plot(rr_transport(prange),grid.zp(prange),'r')
title('transport')
set(gca,'fontsize',fs)
subplot(1,3,2)
plot(rr_relabel(prange),grid.zp(prange),'r')
title('relabel')
set(gca,'fontsize',fs)
subplot(1,3,3)
plot(rr_transport(prange)+rr_relabel(prange),grid.zp(prange),'r')
title('transport + relab')
set(gca,'fontsize',fs)
%pause
%figure(1)



