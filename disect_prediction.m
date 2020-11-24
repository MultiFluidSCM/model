% Disect predicted w2 residual into contributions from different terms

ix =   2:9:9*nz+2;
ixb = 11:9:9*nz+2;
ixt =  2:9:9*nz-7;
ixc = 11:9:9*nz-7;
ikb = 2:nzp;
ikt = 1:nz;
ikc = 2:nz;

% Tendency term
qq = zeros(19,9*nz+6);
if ~switches.c
    qq(10,ix) = qq(10,ix) + m2bar;
    qq(7 ,ix) = qq(7 ,ix) + (w2 - work.w2ubar).*belowr;
    qq(16,ix) = qq(16,ix) + (w2 - work.w2ubar).*abover;
else
    qq(10,2     ) = m2bar(1);
    qq(10,9*nz+2) = m2bar(nzp);
end
rr = - Ndiagmult(qq,xx);
%rr_w1   = rr(1:9:9*nz+1);
rr_tendency   = rr(2:9:9*nz+2);
%rr_eta1 = rr(3:9:9*nz+3);
%rr_eta2 = rr(4:9:9*nz+4);
%rr_q1   = rr(5:9:9*nz+5);
%rr_q2   = rr(6:9:9*nz+6);
%rr_m1   = rr(7:9:9*nz-2);
%rr_m2   = rr(8:9:9*nz-1);
%rr_s    = rr(9:9:9*nz);

% Pressure gradient term + buoyancy
qq = zeros(19,9*nz+6);
qq(8 ,ixc) = qq(8 ,ixc) - adt*m2bar(ikc)./(eos.rhow2(ikc).*dzw(ikc));
qq(17,ixc) = qq(17,ixc) + adt*m2bar(ikc)./(eos.rhow2(ikc).*dzw(ikc));
qq(8 ,ixc) = qq(8 ,ixc) - adt*m2bar(ikc).*dpdz(ikc).*eos.drdpbar2(ikc).*beloww(ikc);
qq(17,ixc) = qq(17,ixc) - adt*m2bar(ikc).*dpdz(ikc).*eos.drdpbar2(ikc).*abovew(ikc);
qq(12,ixc) = qq(12,ixc) - adt*m2bar(ikc).*dpdz(ikc).*eos.drdeta2(ikc);
qq(14,ixc) = qq(14,ixc) - adt*m2bar(ikc).*dpdz(ikc).*eos.drdq2(ikc);
rr = - Ndiagmult(qq,xx);
rr_pgterm   = rr(2:9:9*nz+2);

% Diffusion term
qq = zeros(19,9*nz+6);
qq(1 ,ixc) = qq(1 ,ixc) - adt*work.dDw2dwb(1:nz-1)./dzw(ikc);
qq(10,ixc) = qq(10,ixc) - adt*work.dDw2dwa(1:nz-1)./dzw(ikc);
qq(10,ixc) = qq(10,ixc) + adt*work.dDw2dwb(2:nz  )./dzw(ikc);
qq(19,ixc) = qq(19,ixc) + adt*work.dDw2dwa(2:nz  )./dzw(ikc);
qq(7 ,ixc) = qq(7, ixc) - adt*work.dDw2dm(1:nz-1)./dzw(ikc);
qq(16,ixc) = qq(16,ixc) + adt*work.dDw2dm(2:nz  )./dzw(ikc);
rr = - Ndiagmult(qq,xx);
rr_diffuse   = rr(2:9:9*nz+2);

% Transport term
qq = zeros(19,9*nz+6);
qq(1 ,ixc) = qq(1 ,ixc) - adt*work.dFw2dwb(1:nz-1)./dzw(ikc);
qq(10,ixc) = qq(10,ixc) - adt*work.dFw2dwa(1:nz-1)./dzw(ikc);
qq(10,ixc) = qq(10,ixc) + adt*work.dFw2dwb(2:nz  )./dzw(ikc);
qq(19,ixc) = qq(19,ixc) + adt*work.dFw2dwa(2:nz  )./dzw(ikc);
qq(10,ixc) = qq(10,ixc) + adt*(abovew(2:nz).*belowp(2:nz) + beloww(2:nz).*abovep(1:nz-1)) ...
                            .*m2bar(2:nz).*work.dw2udz(2:nz);
rr = - Ndiagmult(qq,xx);
rr_transport  = rr(2:9:9*nz+2);

% Drag term
qq = zeros(19,9*nz+6);
qq(9 ,ixc) = qq(9 ,ixc) + adt*work.ddragdw1(ikc);
qq(10,ixc) = qq(10,ixc) + adt*work.ddragdw2(ikc);
rr = - Ndiagmult(qq,xx);
rr_drag    = rr(2:9:9*nz+2);

% Relabelling terms
qq = zeros(19,9*nz+6);
qq(10,ixc) = qq(10,ixc) - adt*M21bar(ikc).*relabel.dwhat21dw2(ikc) + adt*M12bar(ikc).*relabel.dwhat12dw2(ikc);
qq(9 ,ixc) = qq(9 ,ixc) - adt*M21bar(ikc).*relabel.dwhat21dw1(ikc) + adt*M12bar(ikc).*relabel.dwhat12dw1(ikc);
qq(6 ,ixb) = qq(6 ,ixb) - adt*belowr(ikb).*(relabel.dM21dm1  .*what21_w2(ikb) - relabel.dM12dm1  .*what12_w2(ikb));
qq(15,ixt) = qq(15,ixt) - adt*abover(ikt).*(relabel.dM21dm1  .*what21_w2(ikt) - relabel.dM12dm1  .*what12_w2(ikt));
qq(7 ,ixb) = qq(7 ,ixb) - adt*belowr(ikb).*(relabel.dM21dm2  .*what21_w2(ikb) - relabel.dM12dm2  .*what12_w2(ikb));
qq(16,ixt) = qq(16,ixt) - adt*abover(ikt).*(relabel.dM21dm2  .*what21_w2(ikt) - relabel.dM12dm2  .*what12_w2(ikt));
qq(9 ,ixb) = qq(9 ,ixb) - adt*(relabel.dM21dw1  .*what21_w2(ikb) - relabel.dM12dw1  .*what12_w2(ikb)).*brap;
qq(9 ,ixt) = qq(9 ,ixt) - adt*(relabel.dM21dw1  .*what21_w2(ikt) - relabel.dM12dw1  .*what12_w2(ikt)).*arbp;
qq(10,ixb) = qq(10,ixb) - adt*(relabel.dM21dw2  .*what21_w2(ikb) - relabel.dM12dw2  .*what12_w2(ikb)).*brap;
qq(10,ixt) = qq(10,ixt) - adt*(relabel.dM21dw2  .*what21_w2(ikt) - relabel.dM12dw2  .*what12_w2(ikt)).*arbp;
qq(11,ixb) = qq(11,ixb) - adt*(relabel.dM21deta1.*what21_w2(ikb) - relabel.dM12deta1.*what12_w2(ikb)).*brap;
qq(11,ixt) = qq(11,ixt) - adt*(relabel.dM21deta1.*what21_w2(ikt) - relabel.dM12deta1.*what12_w2(ikt)).*arbp;
qq(12,ixb) = qq(12,ixb) - adt*(relabel.dM21deta2.*what21_w2(ikb) - relabel.dM12deta2.*what12_w2(ikb)).*brap;
qq(12,ixt) = qq(12,ixt) - adt*(relabel.dM21deta2.*what21_w2(ikt) - relabel.dM12deta2.*what12_w2(ikt)).*arbp;
qq(13,ixb) = qq(13,ixb) - adt*(relabel.dM21dq1  .*what21_w2(ikb) - relabel.dM12dq1  .*what12_w2(ikb)).*brap;
qq(13,ixt) = qq(13,ixt) - adt*(relabel.dM21dq1  .*what21_w2(ikt) - relabel.dM12dq1  .*what12_w2(ikt)).*arbp;
qq(14,ixb) = qq(14,ixb) - adt*(relabel.dM21dq2  .*what21_w2(ikb) - relabel.dM12dq2  .*what12_w2(ikb)).*brap;
qq(14,ixt) = qq(14,ixt) - adt*(relabel.dM21dq2  .*what21_w2(ikt) - relabel.dM12dq2  .*what12_w2(ikt)).*arbp;
rr = - Ndiagmult(qq,xx);
rr_relabel  = rr(2:9:9*nz+2);

rr = rr_transport + rr_diffuse + rr_drag + rr_pgterm + rr_relabel; % + rr_mbar + rr_td;


figure(14)
subplot(2,3,1)
plot(rr_transport(prange),grid.zw(prange),'r')
title('transport')
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
plot(rr_pgterm(prange),grid.zw(prange),'r')
title('pgterm')
set(gca,'fontsize',fs)
subplot(2,3,3)
plot(rr_relabel(prange),grid.zw(prange),'r')
title('relabel')
set(gca,'fontsize',fs)
subplot(2,3,5)
plot(rr_transport(prange)+rr_relabel(prange),grid.zw(prange),'r')
title('trans+relab')
set(gca,'fontsize',fs)
subplot(2,3,6)
plot(rr(prange),grid.zw(prange),'r')
title('total')
set(gca,'fontsize',fs)




ix =   8:9:9*nz-1;
ixb = 17:9:9*nz-1;
ixt =  8:9:9*nz-10;

% Transport term
qq = zeros(19,9*nz+6);
qq(1 ,ix) = qq(1 ,ix) - adt*work.dF2dmb(1:nz)./dzp(1:nz);
qq(10,ix) = qq(10,ix) - adt*work.dF2dma(1:nz)./dzp(1:nz);
qq(10,ix) = qq(10,ix) + adt*work.dF2dmb(2:nzp)./dzp(1:nz);
qq(19,ix) = qq(19,ix) + adt*work.dF2dma(2:nzp)./dzp(1:nz);
qq(4 ,ix) = qq(4 ,ix) - adt*work.dF2dw(1:nz )./dzp;
qq(13,ix) = qq(13,ix) + adt*work.dF2dw(2:nzp)./dzp;
rr = - Ndiagmult(qq,xx);
rr_transport  = rr(8:9:9*nz-1);

% Relabelling terms
qq = zeros(19,9*nz+6);
qq(3 ,ix) = qq(3 ,ix) + adt*(relabel.dM12dw1   - relabel.dM21dw1  ).*belowp;
qq(4 ,ix) = qq(4 ,ix) + adt*(relabel.dM12dw2   - relabel.dM21dw2  ).*belowp;
qq(5 ,ix) = qq(5 ,ix) + adt*(relabel.dM12deta1 - relabel.dM21deta1).*belowp;
qq(6 ,ix) = qq(6 ,ix) + adt*(relabel.dM12deta2 - relabel.dM21deta2).*belowp;
qq(7 ,ix) = qq(7 ,ix) + adt*(relabel.dM12dq1   - relabel.dM21dq1  ).*belowp;
qq(8 ,ix) = qq(8 ,ix) + adt*(relabel.dM12dq2   - relabel.dM21dq2  ).*belowp;
qq(9 ,ix) = qq(9 ,ix) + adt*(relabel.dM12dm1   - relabel.dM21dm1  );
qq(10,ix) = qq(10,ix) + adt*(relabel.dM12dm2   - relabel.dM21dm2  );
qq(12,ix) = qq(12,ix) + adt*(relabel.dM12dw1   - relabel.dM21dw1  ).*abovep;
qq(13,ix) = qq(13,ix) + adt*(relabel.dM12dw2   - relabel.dM21dw2  ).*abovep;
qq(14,ix) = qq(14,ix) + adt*(relabel.dM12deta1 - relabel.dM21deta1).*abovep;
qq(15,ix) = qq(15,ix) + adt*(relabel.dM12deta2 - relabel.dM21deta2).*abovep;
qq(16,ix) = qq(16,ix) + adt*(relabel.dM12dq1   - relabel.dM21dq1  ).*abovep;
qq(17,ix) = qq(17,ix) + adt*(relabel.dM12dq2   - relabel.dM21dq2  ).*abovep;
rr = - Ndiagmult(qq,xx);
rr_relabel  = rr(8:9:9*nz-1);


% Relabelling terms
qq = zeros(19,9*nz+6);
qq(9 ,ix) = qq(9 ,ix) + adt*(relabel.dM12dm1   - relabel.dM21dm1  );
qq(10,ix) = qq(10,ix) + adt*(relabel.dM12dm2   - relabel.dM21dm2  );
rr = - Ndiagmult(qq,xx);
rr_check  = rr(8:9:9*nz-1);


figure(15)
subplot(1,3,1)
plot(rr_transport(prange),grid.zp(prange),'r')
title('pred trans')
set(gca,'fontsize',fs)
subplot(1,3,2)
plot(rr_relabel(prange),grid.zp(prange),'r',rr_check(prange),grid.zp(prange),'b:')
title('pred relab')
set(gca,'fontsize',fs)
subplot(1,3,3)
plot(rr_transport(prange)+rr_relabel(prange),grid.zp(prange),'r')
title('transport + relab')
set(gca,'fontsize',fs)
%pause
%figure(1)







