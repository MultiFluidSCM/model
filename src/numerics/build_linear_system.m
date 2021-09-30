% Build the matrix for w-eta-q-m-p system

% *** The linearization of the entrainment effect of diffusion
% is not yet included ***


% Initialize
cc = zeros(19,9*nz+6);

% Safety parameter in case mbar changes from close to zero to a finite
% value, invalidating the linearization of the w, eta and q equations
% mbar_safety = 1.0e-3;

% Allow for possible relatively large changes in mbar during increments, which
% could invalidate the linearization
%massterm1 = max(m1bar,m1bar + res1mbar); % + mbar_safety;
%massterm2 = max(m2bar,m2bar + res2mbar); % + mbar_safety;

massterm1 = m1bar;
massterm2 = m2bar;
for k = 2:nzp
    massterm1(k) = max(massterm1(k),0.1*massterm1(k-1));
    massterm2(k) = max(massterm2(k),0.1*massterm2(k-1));
end

% w equations

ix =   1:9:9*nz+1;
ixb = 10:9:9*nz+1;
ixt =  1:9:9*nz-8;
ixc = 10:9:9*nz-8;
ikb = 2:nzp;
ikt = 1:nz;
ikc = 2:nz;
% Tendency term
cc(10,ix) = cc(10,ix) + massterm1;
cc(7 ,ix) = cc(7 ,ix) + (w1 - work.w1ubar).*belowr;
cc(16,ix) = cc(16,ix) + (w1 - work.w1ubar).*abover;
if switches.c
    cc(11,ix) = cc(11,ix) + massterm2;
    cc(8 ,ix) = cc(8 ,ix) + (w2 - work.w1ubar).*belowr;
    cc(17,ix) = cc(17,ix) + (w2 - work.w1ubar).*abover;
end
% Pressure gradient term
cc(9 ,ixc) = cc(9 ,ixc) - adt*m1bar(ikc)./(eos.rhow1(ikc).*dzw(ikc));
cc(18,ixc) = cc(18,ixc) + adt*m1bar(ikc)./(eos.rhow1(ikc).*dzw(ikc));
cc(7 ,ixc) = cc(7 ,ixc) + adt*work.nhpg1(ikc).*belowr(ikc);
cc(16,ixc) = cc(16,ixc) + adt*work.nhpg1(ikc).*abover(ikc);
% Buoyancy term
cc(9 ,ixc) = cc(9 ,ixc) - adt*m1bar(ikc).*dpdz(ikc).*eos.drdpbar1(ikc).*beloww(ikc);
cc(18,ixc) = cc(18,ixc) - adt*m1bar(ikc).*dpdz(ikc).*eos.drdpbar1(ikc).*abovew(ikc);
cc(12,ixc) = cc(12,ixc) - adt*m1bar(ikc).*dpdz(ikc).*eos.drdeta1(ikc);
cc(14,ixc) = cc(14,ixc) - adt*m1bar(ikc).*dpdz(ikc).*eos.drdq1(ikc);
% Diffusion term
cc(1 ,ixc) = cc(1 ,ixc) - adt*work.dDw1dwb(1:nz-1)./dzw(ikc);
cc(10,ixc) = cc(10,ixc) - adt*work.dDw1dwa(1:nz-1)./dzw(ikc);
cc(10,ixc) = cc(10,ixc) + adt*work.dDw1dwb(2:nz  )./dzw(ikc);
cc(19,ixc) = cc(19,ixc) + adt*work.dDw1dwa(2:nz  )./dzw(ikc);
cc(7 ,ixc) = cc(7, ixc) - adt*work.dDw1dm(1:nz-1)./dzw(ikc);
cc(16,ixc) = cc(16,ixc) + adt*work.dDw1dm(2:nz  )./dzw(ikc);
% Transport term
cc(1 ,ixc) = cc(1 ,ixc) - adt*work.dFw1dwb(1:nz-1)./dzw(ikc);
cc(10,ixc) = cc(10,ixc) - adt*work.dFw1dwa(1:nz-1)./dzw(ikc);
cc(10,ixc) = cc(10,ixc) + adt*work.dFw1dwb(2:nz  )./dzw(ikc);
cc(19,ixc) = cc(19,ixc) + adt*work.dFw1dwa(2:nz  )./dzw(ikc);
cc(10,ixc) = cc(10,ixc) + adt*(abovew(2:nz).*belowp(2:nz) + beloww(2:nz).*abovep(1:nz-1)) ...
                            .*m1bar(2:nz).*work.dw1udz(2:nz);
% Drag term
cc(10,ixc) = cc(10,ixc) - adt*work.ddragdw1(ikc);
cc(11,ixc) = cc(11,ixc) - adt*work.ddragdw2(ikc);
% Relabelling terms
% ** Only set ixc entries ? **
cc(10,ixc) = cc(10,ixc) - adt*M12bar(ikc).*relabel.dwhat12dw1(ikc) + adt*M21bar(ikc).*relabel.dwhat21dw1(ikc);
cc(11,ixc) = cc(11,ixc) - adt*M12bar(ikc).*relabel.dwhat12dw2(ikc) + adt*M21bar(ikc).*relabel.dwhat21dw2(ikc);
cc(7 ,ixb) = cc(7 ,ixb) - adt*belowr(ikb).*(relabel.dM12dm1  .*what12_w1(ikb) - relabel.dM21dm1  .*what21_w1(ikb));
cc(16,ixt) = cc(16,ixt) - adt*abover(ikt).*(relabel.dM12dm1  .*what12_w1(ikt) - relabel.dM21dm1  .*what21_w1(ikt));
cc(8 ,ixb) = cc(8 ,ixb) - adt*belowr(ikb).*(relabel.dM12dm2  .*what12_w1(ikb) - relabel.dM21dm2  .*what21_w1(ikb));
cc(17,ixt) = cc(17,ixt) - adt*abover(ikt).*(relabel.dM12dm2  .*what12_w1(ikt) - relabel.dM21dm2  .*what21_w1(ikt));
cc(10,ixb) = cc(10,ixb) - adt*(relabel.dM12dw1  .*what12_w1(ikb) - relabel.dM21dw1  .*what21_w1(ikb)).*brap;
cc(10,ixt) = cc(10,ixt) - adt*(relabel.dM12dw1  .*what12_w1(ikt) - relabel.dM21dw1  .*what21_w1(ikt)).*arbp;
cc(11,ixb) = cc(11,ixb) - adt*(relabel.dM12dw2  .*what12_w1(ikb) - relabel.dM21dw2  .*what21_w1(ikb)).*brap;
cc(11,ixt) = cc(11,ixt) - adt*(relabel.dM12dw2  .*what12_w1(ikt) - relabel.dM21dw2  .*what21_w1(ikt)).*arbp;
cc(12,ixb) = cc(12,ixb) - adt*(relabel.dM12deta1.*what12_w1(ikb) - relabel.dM21deta1.*what21_w1(ikb)).*brap;
cc(12,ixt) = cc(12,ixt) - adt*(relabel.dM12deta1.*what12_w1(ikt) - relabel.dM21deta1.*what21_w1(ikt)).*arbp;
cc(13,ixb) = cc(13,ixb) - adt*(relabel.dM12deta2.*what12_w1(ikb) - relabel.dM21deta2.*what21_w1(ikb)).*brap;
cc(13,ixt) = cc(13,ixt) - adt*(relabel.dM12deta2.*what12_w1(ikt) - relabel.dM21deta2.*what21_w1(ikt)).*arbp;
cc(14,ixb) = cc(14,ixb) - adt*(relabel.dM12dq1  .*what12_w1(ikb) - relabel.dM21dq1  .*what21_w1(ikb)).*brap;
cc(14,ixt) = cc(14,ixt) - adt*(relabel.dM12dq1  .*what12_w1(ikt) - relabel.dM21dq1  .*what21_w1(ikt)).*arbp;
cc(15,ixb) = cc(15,ixb) - adt*(relabel.dM12dq2  .*what12_w1(ikb) - relabel.dM21dq2  .*what21_w1(ikb)).*brap;
cc(15,ixt) = cc(15,ixt) - adt*(relabel.dM12dq2  .*what12_w1(ikt) - relabel.dM21dq2  .*what21_w1(ikt)).*arbp;

ix =   2:9:9*nz+2;
ixb = 11:9:9*nz+2;
ixt =  2:9:9*nz-7;
ixc = 11:9:9*nz-7;
ikb = 2:nzp;
ikt = 1:nz;
ikc = 2:nz;
% Tendency term
if ~switches.c
    cc(10,ix) = cc(10,ix) + massterm2;
    cc(7 ,ix) = cc(7 ,ix) + (w2 - work.w2ubar).*belowr;
    cc(16,ix) = cc(16,ix) + (w2 - work.w2ubar).*abover;
else
    cc(10,2     ) = massterm2(1);
    cc(10,9*nz+2) = massterm2(nzp);
end
% Pressure gradient term
cc(8 ,ixc) = cc(8 ,ixc) - adt*m2bar(ikc)./(eos.rhow2(ikc).*dzw(ikc));
cc(17,ixc) = cc(17,ixc) + adt*m2bar(ikc)./(eos.rhow2(ikc).*dzw(ikc));
cc(7 ,ixc) = cc(7 ,ixc) + adt*work.nhpg2(ikc).*beloww(ikc);
cc(16,ixc) = cc(16,ixc) + adt*work.nhpg2(ikc).*abovew(ikc);
% Buoyancy term
cc(8 ,ixc) = cc(8 ,ixc) - adt*m2bar(ikc).*dpdz(ikc).*eos.drdpbar2(ikc).*beloww(ikc);
cc(17,ixc) = cc(17,ixc) - adt*m2bar(ikc).*dpdz(ikc).*eos.drdpbar2(ikc).*abovew(ikc);
cc(12,ixc) = cc(12,ixc) - adt*m2bar(ikc).*dpdz(ikc).*eos.drdeta2(ikc);
cc(14,ixc) = cc(14,ixc) - adt*m2bar(ikc).*dpdz(ikc).*eos.drdq2(ikc);
% Diffusion term
cc(1 ,ixc) = cc(1 ,ixc) - adt*work.dDw2dwb(1:nz-1)./dzw(ikc);
cc(10,ixc) = cc(10,ixc) - adt*work.dDw2dwa(1:nz-1)./dzw(ikc);
cc(10,ixc) = cc(10,ixc) + adt*work.dDw2dwb(2:nz  )./dzw(ikc);
cc(19,ixc) = cc(19,ixc) + adt*work.dDw2dwa(2:nz  )./dzw(ikc);
cc(7 ,ixc) = cc(7, ixc) - adt*work.dDw2dm(1:nz-1)./dzw(ikc);
cc(16,ixc) = cc(16,ixc) + adt*work.dDw2dm(2:nz  )./dzw(ikc);
% Transport term
cc(1 ,ixc) = cc(1 ,ixc) - adt*work.dFw2dwb(1:nz-1)./dzw(ikc);
cc(10,ixc) = cc(10,ixc) - adt*work.dFw2dwa(1:nz-1)./dzw(ikc);
cc(10,ixc) = cc(10,ixc) + adt*work.dFw2dwb(2:nz  )./dzw(ikc);
cc(19,ixc) = cc(19,ixc) + adt*work.dFw2dwa(2:nz  )./dzw(ikc);
cc(10,ixc) = cc(10,ixc) + adt*(abovew(2:nz).*belowp(2:nz) + beloww(2:nz).*abovep(1:nz-1)) ...
                            .*m2bar(2:nz).*work.dw2udz(2:nz);
% Drag term
cc(9 ,ixc) = cc(9 ,ixc) + adt*work.ddragdw1(ikc);
cc(10,ixc) = cc(10,ixc) + adt*work.ddragdw2(ikc);
% Relabelling terms
cc(10,ixc) = cc(10,ixc) - adt*M21bar(ikc).*relabel.dwhat21dw2(ikc) + adt*M12bar(ikc).*relabel.dwhat12dw2(ikc);
cc(9 ,ixc) = cc(9 ,ixc) - adt*M21bar(ikc).*relabel.dwhat21dw1(ikc) + adt*M12bar(ikc).*relabel.dwhat12dw1(ikc);
cc(6 ,ixb) = cc(6 ,ixb) - adt*belowr(ikb).*(relabel.dM21dm1  .*what21_w2(ikb) - relabel.dM12dm1  .*what12_w2(ikb));
cc(15,ixt) = cc(15,ixt) - adt*abover(ikt).*(relabel.dM21dm1  .*what21_w2(ikt) - relabel.dM12dm1  .*what12_w2(ikt));
cc(7 ,ixb) = cc(7 ,ixb) - adt*belowr(ikb).*(relabel.dM21dm2  .*what21_w2(ikb) - relabel.dM12dm2  .*what12_w2(ikb));
cc(16,ixt) = cc(16,ixt) - adt*abover(ikt).*(relabel.dM21dm2  .*what21_w2(ikt) - relabel.dM12dm2  .*what12_w2(ikt));
cc(9 ,ixb) = cc(9 ,ixb) - adt*(relabel.dM21dw1  .*what21_w2(ikb) - relabel.dM12dw1  .*what12_w2(ikb)).*brap;
cc(9 ,ixt) = cc(9 ,ixt) - adt*(relabel.dM21dw1  .*what21_w2(ikt) - relabel.dM12dw1  .*what12_w2(ikt)).*arbp;
cc(10,ixb) = cc(10,ixb) - adt*(relabel.dM21dw2  .*what21_w2(ikb) - relabel.dM12dw2  .*what12_w2(ikb)).*brap;
cc(10,ixt) = cc(10,ixt) - adt*(relabel.dM21dw2  .*what21_w2(ikt) - relabel.dM12dw2  .*what12_w2(ikt)).*arbp;
cc(11,ixb) = cc(11,ixb) - adt*(relabel.dM21deta1.*what21_w2(ikb) - relabel.dM12deta1.*what12_w2(ikb)).*brap;
cc(11,ixt) = cc(11,ixt) - adt*(relabel.dM21deta1.*what21_w2(ikt) - relabel.dM12deta1.*what12_w2(ikt)).*arbp;
cc(12,ixb) = cc(12,ixb) - adt*(relabel.dM21deta2.*what21_w2(ikb) - relabel.dM12deta2.*what12_w2(ikb)).*brap;
cc(12,ixt) = cc(12,ixt) - adt*(relabel.dM21deta2.*what21_w2(ikt) - relabel.dM12deta2.*what12_w2(ikt)).*arbp;
cc(13,ixb) = cc(13,ixb) - adt*(relabel.dM21dq1  .*what21_w2(ikb) - relabel.dM12dq1  .*what12_w2(ikb)).*brap;
cc(13,ixt) = cc(13,ixt) - adt*(relabel.dM21dq1  .*what21_w2(ikt) - relabel.dM12dq1  .*what12_w2(ikt)).*arbp;
cc(14,ixb) = cc(14,ixb) - adt*(relabel.dM21dq2  .*what21_w2(ikb) - relabel.dM12dq2  .*what12_w2(ikb)).*brap;
cc(14,ixt) = cc(14,ixt) - adt*(relabel.dM21dq2  .*what21_w2(ikt) - relabel.dM12dq2  .*what12_w2(ikt)).*arbp;

% eta equations

ix =   3:9:9*nz+3;
ixb = 12:9:9*nz+3;
ixt =  3:9:9*nz-6;
ixm = 12:9:9*nz-6;
ikb = 2:nzp;
ikt = 1:nz;
ikc = 2:nz;
% Tendency term
cc(10,ix) = cc(10,ix) + massterm1;
cc(5 ,ix) = cc(5 ,ix) + (eta1 - work.eta1ubar).*belowr;
cc(14,ix) = cc(14,ix) + (eta1 - work.eta1ubar).*abover;
if switches.c
    cc(11,ix) = cc(11,ix) + massterm2;
    cc(6 ,ix) = cc(6 ,ix) + (eta2 - work.eta1ubar).*belowr;
    cc(15,ix) = cc(15,ix) + (eta2 - work.eta1ubar).*abover;
end
% Diffusion term
cc(1 ,ixb) = cc(1 ,ixb) - adt*work.dDeta1detab./dzw(ikb);
cc(10,ixb) = cc(10,ixb) - adt*work.dDeta1detaa./dzw(ikb);
cc(10,ixt) = cc(10,ixt) + adt*work.dDeta1detab./dzw(ikt);
cc(19,ixt) = cc(19,ixt) + adt*work.dDeta1detaa./dzw(ikt);
cc(5 ,ixb) = cc(5, ixb) - adt*work.dDeta1dm./dzw(ikb);
cc(14,ixt) = cc(14,ixt) + adt*work.dDeta1dm./dzw(ikt);
% Transport term
cc(1 ,ixb) = cc(1 ,ixb) - adt*work.dFeta1detab./dzw(ikb);
cc(10,ixb) = cc(10,ixb) - adt*work.dFeta1detaa./dzw(ikb);
cc(10,ixt) = cc(10,ixt) + adt*work.dFeta1detab./dzw(ikt);
cc(19,ixt) = cc(19,ixt) + adt*work.dFeta1detaa./dzw(ikt);
cc(8 ,ixm) = cc(8 ,ixm) + adt*(abovew(2:nz).*belowp(2:nz) + beloww(2:nz).*abovep(1:nz-1)) ...
                            .*m1bar(2:nz).*work.deta1udz(2:nz);
% Relabelling terms
cc(10,ix) = cc(10,ix) - adt*M12bar.*relabel.detahat12deta1 + adt*M21bar.*relabel.detahat21deta1;
cc(11,ix) = cc(11,ix) - adt*M12bar.*relabel.detahat12deta2 + adt*M21bar.*relabel.detahat21deta2;
cc(5 ,ixb) = cc(5 ,ixb) - adt*belowr(ikb).*(relabel.dM12dm1  .*etahat12_eta1(ikb) ...
                                          - relabel.dM21dm1  .*etahat21_eta1(ikb));
cc(14,ixt) = cc(14,ixt) - adt*abover(ikt).*(relabel.dM12dm1  .*etahat12_eta1(ikt) ...
                                          - relabel.dM21dm1  .*etahat21_eta1(ikt));
cc(6 ,ixb) = cc(6 ,ixb) - adt*belowr(ikb).*(relabel.dM12dm2  .*etahat12_eta1(ikb) ...
                                          - relabel.dM21dm2  .*etahat21_eta1(ikb));
cc(15,ixt) = cc(15,ixt) - adt*abover(ikt).*(relabel.dM12dm2  .*etahat12_eta1(ikt) ...
                                          - relabel.dM21dm2  .*etahat21_eta1(ikt));
cc(8 ,ixb) = cc(8 ,ixb) - adt*(relabel.dM12dw1  .*etahat12_eta1(ikb) - relabel.dM21dw1  .*etahat21_eta1(ikb)).*brap;
cc(8 ,ixt) = cc(8 ,ixt) - adt*(relabel.dM12dw1  .*etahat12_eta1(ikt) - relabel.dM21dw1  .*etahat21_eta1(ikt)).*arbp;
cc(9 ,ixb) = cc(9 ,ixb) - adt*(relabel.dM12dw2  .*etahat12_eta1(ikb) - relabel.dM21dw2  .*etahat21_eta1(ikb)).*brap;
cc(9 ,ixt) = cc(9 ,ixt) - adt*(relabel.dM12dw2  .*etahat12_eta1(ikt) - relabel.dM21dw2  .*etahat21_eta1(ikt)).*arbp;
cc(10,ixb) = cc(10,ixb) - adt*(relabel.dM12deta1.*etahat12_eta1(ikb) - relabel.dM21deta1.*etahat21_eta1(ikb)).*brap;
cc(10,ixt) = cc(10,ixt) - adt*(relabel.dM12deta1.*etahat12_eta1(ikt) - relabel.dM21deta1.*etahat21_eta1(ikt)).*arbp;
cc(11,ixb) = cc(11,ixb) - adt*(relabel.dM12deta2.*etahat12_eta1(ikb) - relabel.dM21deta2.*etahat21_eta1(ikb)).*brap;
cc(11,ixt) = cc(11,ixt) - adt*(relabel.dM12deta2.*etahat12_eta1(ikt) - relabel.dM21deta2.*etahat21_eta1(ikt)).*arbp;
cc(12,ixb) = cc(12,ixb) - adt*(relabel.dM12dq1  .*etahat12_eta1(ikb) - relabel.dM21dq1  .*etahat21_eta1(ikb)).*brap;
cc(12,ixt) = cc(12,ixt) - adt*(relabel.dM12dq1  .*etahat12_eta1(ikt) - relabel.dM21dq1  .*etahat21_eta1(ikt)).*arbp;
cc(13,ixb) = cc(13,ixb) - adt*(relabel.dM12dq2  .*etahat12_eta1(ikb) - relabel.dM21dq2  .*etahat21_eta1(ikb)).*brap;
cc(13,ixt) = cc(13,ixt) - adt*(relabel.dM12dq2  .*etahat12_eta1(ikt) - relabel.dM21dq2  .*etahat21_eta1(ikt)).*arbp;

ix =   4:9:9*nz+4;
ixb = 13:9:9*nz+4;
ixt =  4:9:9*nz-5;
ixm = 13:9:9*nz-5;
ikb = 2:nzp;
ikt = 1:nz;
ikc = 2:nz;
% Tendency term
if ~switches.c
    cc(10,ix) = cc(10,ix) + massterm2;
    cc(5 ,ix) = cc(5 ,ix) + (eta2 - work.eta2ubar).*belowr;
    cc(14,ix) = cc(14,ix) + (eta2 - work.eta2ubar).*abover;
else
    cc(10,4 ) = massterm2(1);
end
% Diffusion term
cc(1 ,ixb) = cc(1 ,ixb) - adt*work.dDeta2detab./dzw(ikb);
cc(10,ixb) = cc(10,ixb) - adt*work.dDeta2detaa./dzw(ikb);
cc(10,ixt) = cc(10,ixt) + adt*work.dDeta2detab./dzw(ikt);
cc(19,ixt) = cc(19,ixt) + adt*work.dDeta2detaa./dzw(ikt);
cc(5 ,ixb) = cc(5, ixb) - adt*work.dDeta2dm./dzw(ikb);
cc(14,ixt) = cc(14,ixt) + adt*work.dDeta2dm./dzw(ikt);
% Transport term
cc(1 ,ixb) = cc(1 ,ixb) - adt*work.dFeta2detab./dzw(ikb);
cc(10,ixb) = cc(10,ixb) - adt*work.dFeta2detaa./dzw(ikb);
cc(10,ixt) = cc(10,ixt) + adt*work.dFeta2detab./dzw(ikt);
cc(19,ixt) = cc(19,ixt) + adt*work.dFeta2detaa./dzw(ikt);
cc(8 ,ixm) = cc(8 ,ixm) + adt*(abovew(2:nz).*belowp(2:nz) + beloww(2:nz).*abovep(1:nz-1)) ...
                            .*m2bar(2:nz).*work.deta2udz(2:nz);
% Relabelling terms
cc(10,ix) = cc(10,ix) - adt*M21bar.*relabel.detahat21deta2 + adt*M12bar.*relabel.detahat12deta2;
cc(9 ,ix) = cc(9 ,ix) - adt*M21bar.*relabel.detahat21deta1 + adt*M12bar.*relabel.detahat12deta1;
cc(4 ,ixb) = cc(4 ,ixb) - adt*belowr(ikb).*(relabel.dM21dm1  .*etahat21_eta2(ikb) - relabel.dM12dm1  .*etahat12_eta2(ikb));
cc(13,ixt) = cc(13,ixt) - adt*abover(ikt).*(relabel.dM21dm1  .*etahat21_eta2(ikt) - relabel.dM12dm1  .*etahat12_eta2(ikt));
cc(5 ,ixb) = cc(5 ,ixb) - adt*belowr(ikb).*(relabel.dM21dm2  .*etahat21_eta2(ikb) - relabel.dM12dm2  .*etahat12_eta2(ikb));
cc(14,ixt) = cc(14,ixt) - adt*abover(ikt).*(relabel.dM21dm2  .*etahat21_eta2(ikt) - relabel.dM12dm2  .*etahat12_eta2(ikt));
cc(7 ,ixb) = cc(7 ,ixb) - adt*(relabel.dM21dw1  .*etahat21_eta2(ikb) - relabel.dM12dw1  .*etahat12_eta2(ikb)).*brap;
cc(7 ,ixt) = cc(7 ,ixt) - adt*(relabel.dM21dw1  .*etahat21_eta2(ikt) - relabel.dM12dw1  .*etahat12_eta2(ikt)).*arbp;
cc(8 ,ixb) = cc(8 ,ixb) - adt*(relabel.dM21dw2  .*etahat21_eta2(ikb) - relabel.dM12dw2  .*etahat12_eta2(ikb)).*brap;
cc(8 ,ixt) = cc(8 ,ixt) - adt*(relabel.dM21dw2  .*etahat21_eta2(ikt) - relabel.dM12dw2  .*etahat12_eta2(ikt)).*arbp;
cc(9 ,ixb) = cc(9 ,ixb) - adt*(relabel.dM21deta1.*etahat21_eta2(ikb) - relabel.dM12deta1.*etahat12_eta2(ikb)).*brap;
cc(9 ,ixt) = cc(9 ,ixt) - adt*(relabel.dM21deta1.*etahat21_eta2(ikt) - relabel.dM12deta1.*etahat12_eta2(ikt)).*arbp;
cc(10,ixb) = cc(10,ixb) - adt*(relabel.dM21deta2.*etahat21_eta2(ikb) - relabel.dM12deta2.*etahat12_eta2(ikb)).*brap;
cc(10,ixt) = cc(10,ixt) - adt*(relabel.dM21deta2.*etahat21_eta2(ikt) - relabel.dM12deta2.*etahat12_eta2(ikt)).*arbp;
cc(11,ixb) = cc(11,ixb) - adt*(relabel.dM21dq1  .*etahat21_eta2(ikb) - relabel.dM12dq1  .*etahat12_eta2(ikb)).*brap;
cc(11,ixt) = cc(11,ixt) - adt*(relabel.dM21dq1  .*etahat21_eta2(ikt) - relabel.dM12dq1  .*etahat12_eta2(ikt)).*arbp;
cc(12,ixb) = cc(12,ixb) - adt*(relabel.dM21dq2  .*etahat21_eta2(ikb) - relabel.dM12dq2  .*etahat12_eta2(ikb)).*brap;
cc(12,ixt) = cc(12,ixt) - adt*(relabel.dM21dq2  .*etahat21_eta2(ikt) - relabel.dM12dq2  .*etahat12_eta2(ikt)).*arbp;


% q equations

ix =   5:9:9*nz+5;
ixb = 14:9:9*nz+5;
ixt =  5:9:9*nz-4;
ixm = 14:9:9*nz-4;
% Tendency term
cc(10,ix) = cc(10,ix) + massterm1;
cc(3 ,ix) = cc(3 ,ix) + (q1 - work.q1ubar).*belowr;
cc(12,ix) = cc(12,ix) + (q1 - work.q1ubar).*abover;
if switches.c
    cc(11,ix) = cc(11,ix) + massterm2;
    cc(4 ,ix) = cc(4 ,ix) + (q2 - work.q1ubar).*belowr;
    cc(13,ix) = cc(13,ix) + (q2 - work.q1ubar).*abover;
end
% Diffusion term
cc(1 ,ixb) = cc(1 ,ixb) - adt*work.dDq1dqb./dzw(ikb);
cc(10,ixb) = cc(10,ixb) - adt*work.dDq1dqa./dzw(ikb);
cc(10,ixt) = cc(10,ixt) + adt*work.dDq1dqb./dzw(ikt);
cc(19,ixt) = cc(19,ixt) + adt*work.dDq1dqa./dzw(ikt);
cc(3 ,ixb) = cc(3, ixb) - adt*work.dDq1dm./dzw(ikb);
cc(12,ixt) = cc(12,ixt) + adt*work.dDq1dm./dzw(ikt);
% Transport term
cc(1 ,ixb) = cc(1 ,ixb) - adt*work.dFq1dqb./dzw(ikb);
cc(10,ixb) = cc(10,ixb) - adt*work.dFq1dqa./dzw(ikb);
cc(10,ixt) = cc(10,ixt) + adt*work.dFq1dqb./dzw(ikt);
cc(19,ixt) = cc(19,ixt) + adt*work.dFq1dqa./dzw(ikt);
cc(6 ,ixm) = cc(6 ,ixm) + adt*(abovew(2:nz).*belowp(2:nz) + beloww(2:nz).*abovep(1:nz-1)) ...
                            .*m1bar(2:nz).*work.dq1udz(2:nz);
% Relabelling terms
cc(10,ix) = cc(10,ix) - adt*M12bar.*relabel.dqhat12dq1 + adt*M21bar.*relabel.dqhat21dq1;
cc(11,ix) = cc(11,ix) - adt*M12bar.*relabel.dqhat12dq2 + adt*M21bar.*relabel.dqhat21dq2;
cc(3 ,ixb) = cc(3 ,ixb) - adt*belowr(ikb).*(relabel.dM12dm1  .*qhat12_q1(ikb) ...
                                          - relabel.dM21dm1  .*qhat21_q1(ikb));
cc(12,ixt) = cc(12,ixt) - adt*abover(ikt).*(relabel.dM12dm1  .*qhat12_q1(ikt) ...
                                          - relabel.dM21dm1  .*qhat21_q1(ikt));
cc(4 ,ixb) = cc(4 ,ixb) - adt*belowr(ikb).*(relabel.dM12dm2  .*qhat12_q1(ikb) ...
                                          - relabel.dM21dm2  .*qhat21_q1(ikb));
cc(13,ixt) = cc(13,ixt) - adt*abover(ikt).*(relabel.dM12dm2  .*qhat12_q1(ikt) ...
                                          - relabel.dM21dm2  .*qhat21_q1(ikt));
cc(6 ,ixb) = cc(6 ,ixb) - adt*(relabel.dM12dw1  .*qhat12_q1(ikb) - relabel.dM21dw1  .*qhat21_q1(ikb)).*brap;
cc(6 ,ixt) = cc(6 ,ixt) - adt*(relabel.dM12dw1  .*qhat12_q1(ikt) - relabel.dM21dw1  .*qhat21_q1(ikt)).*arbp;
cc(7 ,ixb) = cc(7 ,ixb) - adt*(relabel.dM12dw2  .*qhat12_q1(ikb) - relabel.dM21dw2  .*qhat21_q1(ikb)).*brap;
cc(7 ,ixt) = cc(7 ,ixt) - adt*(relabel.dM12dw2  .*qhat12_q1(ikt) - relabel.dM21dw2  .*qhat21_q1(ikt)).*arbp;
cc(8 ,ixb) = cc(8 ,ixb) - adt*(relabel.dM12deta1.*qhat12_q1(ikb) - relabel.dM21deta1.*qhat21_q1(ikb)).*brap;
cc(8 ,ixt) = cc(8 ,ixt) - adt*(relabel.dM12deta1.*qhat12_q1(ikt) - relabel.dM21deta1.*qhat21_q1(ikt)).*arbp;
cc(9 ,ixb) = cc(9 ,ixb) - adt*(relabel.dM12deta2.*qhat12_q1(ikb) - relabel.dM21deta2.*qhat21_q1(ikb)).*brap;
cc(9 ,ixt) = cc(9 ,ixt) - adt*(relabel.dM12deta2.*qhat12_q1(ikt) - relabel.dM21deta2.*qhat21_q1(ikt)).*arbp;
cc(10,ixb) = cc(10,ixb) - adt*(relabel.dM12dq1  .*qhat12_q1(ikb) - relabel.dM21dq1  .*qhat21_q1(ikb)).*brap;
cc(10,ixt) = cc(10,ixt) - adt*(relabel.dM12dq1  .*qhat12_q1(ikt) - relabel.dM21dq1  .*qhat21_q1(ikt)).*arbp;
cc(11,ixb) = cc(11,ixb) - adt*(relabel.dM12dq2  .*qhat12_q1(ikb) - relabel.dM21dq2  .*qhat21_q1(ikb)).*brap;
cc(11,ixt) = cc(11,ixt) - adt*(relabel.dM12dq2  .*qhat12_q1(ikt) - relabel.dM21dq2  .*qhat21_q1(ikt)).*arbp;

ix =   6:9:9*nz+6;
ixb = 15:9:9*nz+6;
ixt =  6:9:9*nz-3;
ixm = 15:9:9*nz-3;
% Tendency term
if ~switches.c
    cc(10,ix) = cc(10,ix) + massterm2;
    cc(3 ,ix) = cc(3 ,ix) + (q2 - work.q2ubar).*belowr;
    cc(12,ix) = cc(12,ix) + (q2 - work.q2ubar).*abover;
else
    cc(10,6 ) = massterm2(1);
end
% Diffusion term
cc(1 ,ixb) = cc(1 ,ixb) - adt*work.dDq2dqb./dzw(ikb);
cc(10,ixb) = cc(10,ixb) - adt*work.dDq2dqa./dzw(ikb);
cc(10,ixt) = cc(10,ixt) + adt*work.dDq2dqb./dzw(ikt);
cc(19,ixt) = cc(19,ixt) + adt*work.dDq2dqa./dzw(ikt);
cc(3 ,ixb) = cc(3, ixb) - adt*work.dDq2dm./dzw(ikb);
cc(12,ixt) = cc(12,ixt) + adt*work.dDq2dm./dzw(ikt);
% Transport term
cc(1 ,ixb) = cc(1 ,ixb) - adt*work.dFq2dqb./dzw(2:nzp);
cc(10,ixb) = cc(10,ixb) - adt*work.dFq2dqa./dzw(2:nzp);
cc(10,ixt) = cc(10,ixt) + adt*work.dFq2dqb./dzw(1:nz);
cc(19,ixt) = cc(19,ixt) + adt*work.dFq2dqa./dzw(1:nz);
cc(6 ,ixm) = cc(6 ,ixm) + adt*(abovew(2:nz).*belowp(2:nz) + beloww(2:nz).*abovep(1:nz-1)) ...
                            .*m2bar(2:nz).*work.dq2udz(2:nz);
% Relabelling terms
cc(10,ix) = cc(10,ix) - adt*M21bar.*relabel.dqhat21dq2 + adt*M12bar.*relabel.dqhat12dq2;
cc(9 ,ix) = cc(9 ,ix) - adt*M21bar.*relabel.dqhat21dq1 + adt*M12bar.*relabel.dqhat12dq1;
cc(2 ,ixb) = cc(2 ,ixb) - adt*belowr(ikb).*(relabel.dM21dm1  .*qhat21_q2(ikb) - relabel.dM12dm1  .*qhat12_q2(ikb));
cc(11,ixt) = cc(11,ixt) - adt*abover(ikt).*(relabel.dM21dm1  .*qhat21_q2(ikt) - relabel.dM12dm1  .*qhat12_q2(ikt));
cc(3 ,ixb) = cc(3 ,ixb) - adt*belowr(ikb).*(relabel.dM21dm2  .*qhat21_q2(ikb) - relabel.dM12dm2  .*qhat12_q2(ikb));
cc(12,ixt) = cc(12,ixt) - adt*abover(ikt).*(relabel.dM21dm2  .*qhat21_q2(ikt) - relabel.dM12dm2  .*qhat12_q2(ikt));
cc(5 ,ixb) = cc(5 ,ixb) - adt*(relabel.dM21dw1  .*qhat21_q2(ikb) - relabel.dM12dw1  .*qhat12_q2(ikb)).*brap;
cc(5 ,ixt) = cc(5 ,ixt) - adt*(relabel.dM21dw1  .*qhat21_q2(ikt) - relabel.dM12dw1  .*qhat12_q2(ikt)).*arbp;
cc(6 ,ixb) = cc(6 ,ixb) - adt*(relabel.dM21dw2  .*qhat21_q2(ikb) - relabel.dM12dw2  .*qhat12_q2(ikb)).*brap;
cc(6 ,ixt) = cc(6 ,ixt) - adt*(relabel.dM21dw2  .*qhat21_q2(ikt) - relabel.dM12dw2  .*qhat12_q2(ikt)).*arbp;
cc(7 ,ixb) = cc(7 ,ixb) - adt*(relabel.dM21deta1.*qhat21_q2(ikb) - relabel.dM12deta1.*qhat12_q2(ikb)).*brap;
cc(7 ,ixt) = cc(7 ,ixt) - adt*(relabel.dM21deta1.*qhat21_q2(ikt) - relabel.dM12deta1.*qhat12_q2(ikt)).*arbp;
cc(8 ,ixb) = cc(8 ,ixb) - adt*(relabel.dM21deta2.*qhat21_q2(ikb) - relabel.dM12deta2.*qhat12_q2(ikb)).*brap;
cc(8 ,ixt) = cc(8 ,ixt) - adt*(relabel.dM21deta2.*qhat21_q2(ikt) - relabel.dM12deta2.*qhat12_q2(ikt)).*arbp;
cc(9 ,ixb) = cc(9 ,ixb) - adt*(relabel.dM21dq1  .*qhat21_q2(ikb) - relabel.dM12dq1  .*qhat12_q2(ikb)).*brap;
cc(9 ,ixt) = cc(9 ,ixt) - adt*(relabel.dM21dq1  .*qhat21_q2(ikt) - relabel.dM12dq1  .*qhat12_q2(ikt)).*arbp;
cc(10,ixb) = cc(10,ixb) - adt*(relabel.dM21dq2  .*qhat21_q2(ikb) - relabel.dM12dq2  .*qhat12_q2(ikb)).*brap;
cc(10,ixt) = cc(10,ixt) - adt*(relabel.dM21dq2  .*qhat21_q2(ikt) - relabel.dM12dq2  .*qhat12_q2(ikt)).*arbp;


% Mass equations

ix =   7:9:9*nz-2;
ixb = 16:9:9*nz-2;
ixt =  7:9:9*nz-11;
% Tendency term
cc(10,ix) = 1;
if switches.c
    cc(11,ix) = cc(11,ix) + 1;
end
% Transport term
cc(1 ,ix) = cc(1 ,ix) - adt*work.dF1dmb(1:nz )./dzp(1:nz);
cc(10,ix) = cc(10,ix) - adt*work.dF1dma(1:nz )./dzp(1:nz);
cc(10,ix) = cc(10,ix) + adt*work.dF1dmb(2:nzp)./dzp(1:nz);
cc(19,ix) = cc(19,ix) + adt*work.dF1dma(2:nzp)./dzp(1:nz);
cc(4 ,ix) = cc(4 ,ix) - adt*work.dF1dw(1:nz )./dzp;
cc(13,ix) = cc(13,ix) + adt*work.dF1dw(2:nzp)./dzp;
% Relabelling terms
cc(4 ,ix) = cc(4 ,ix) + adt*(relabel.dM21dw1   - relabel.dM12dw1  ).*belowp;
cc(5 ,ix) = cc(5 ,ix) + adt*(relabel.dM21dw2   - relabel.dM12dw2  ).*belowp;
cc(6 ,ix) = cc(6 ,ix) + adt*(relabel.dM21deta1 - relabel.dM12deta1).*belowp;
cc(7 ,ix) = cc(7 ,ix) + adt*(relabel.dM21deta2 - relabel.dM12deta2).*belowp;
cc(8 ,ix) = cc(8 ,ix) + adt*(relabel.dM21dq1   - relabel.dM12dq1  ).*belowp;
cc(9 ,ix) = cc(9 ,ix) + adt*(relabel.dM21dq2   - relabel.dM12dq2  ).*belowp;
cc(10,ix) = cc(10,ix) + adt*(relabel.dM21dm1   - relabel.dM12dm1  );
cc(11,ix) = cc(11,ix) + adt*(relabel.dM21dm2   - relabel.dM12dm2  );
cc(13,ix) = cc(13,ix) + adt*(relabel.dM21dw1   - relabel.dM12dw1  ).*abovep;
cc(14,ix) = cc(14,ix) + adt*(relabel.dM21dw2   - relabel.dM12dw2  ).*abovep;
cc(15,ix) = cc(15,ix) + adt*(relabel.dM21deta1 - relabel.dM12deta1).*abovep;
cc(16,ix) = cc(16,ix) + adt*(relabel.dM21deta2 - relabel.dM12deta2).*abovep;
cc(17,ix) = cc(17,ix) + adt*(relabel.dM21dq1   - relabel.dM12dq1  ).*abovep;
cc(18,ix) = cc(18,ix) + adt*(relabel.dM21dq2   - relabel.dM12dq2  ).*abovep;

ix =   8:9:9*nz-1;
ixb = 17:9:9*nz-1;
ixt =  8:9:9*nz-10;
% Tendency term
if ~switches.c
    cc(10,ix) = 1;
end
% Transport term
cc(1 ,ix) = cc(1 ,ix) - adt*work.dF2dmb(1:nz)./dzp(1:nz);
cc(10,ix) = cc(10,ix) - adt*work.dF2dma(1:nz)./dzp(1:nz);
cc(10,ix) = cc(10,ix) + adt*work.dF2dmb(2:nzp)./dzp(1:nz);
cc(19,ix) = cc(19,ix) + adt*work.dF2dma(2:nzp)./dzp(1:nz);
cc(4 ,ix) = cc(4 ,ix) - adt*work.dF2dw(1:nz )./dzp;
cc(13,ix) = cc(13,ix) + adt*work.dF2dw(2:nzp)./dzp;
% Relabelling terms
cc(3 ,ix) = cc(3 ,ix) + adt*(relabel.dM12dw1   - relabel.dM21dw1  ).*belowp;
cc(4 ,ix) = cc(4 ,ix) + adt*(relabel.dM12dw2   - relabel.dM21dw2  ).*belowp;
cc(5 ,ix) = cc(5 ,ix) + adt*(relabel.dM12deta1 - relabel.dM21deta1).*belowp;
cc(6 ,ix) = cc(6 ,ix) + adt*(relabel.dM12deta2 - relabel.dM21deta2).*belowp;
cc(7 ,ix) = cc(7 ,ix) + adt*(relabel.dM12dq1   - relabel.dM21dq1  ).*belowp;
cc(8 ,ix) = cc(8 ,ix) + adt*(relabel.dM12dq2   - relabel.dM21dq2  ).*belowp;
cc(9 ,ix) = cc(9 ,ix) + adt*(relabel.dM12dm1   - relabel.dM21dm1  );
cc(10,ix) = cc(10,ix) + adt*(relabel.dM12dm2   - relabel.dM21dm2  );
cc(12,ix) = cc(12,ix) + adt*(relabel.dM12dw1   - relabel.dM21dw1  ).*abovep;
cc(13,ix) = cc(13,ix) + adt*(relabel.dM12dw2   - relabel.dM21dw2  ).*abovep;
cc(14,ix) = cc(14,ix) + adt*(relabel.dM12deta1 - relabel.dM21deta1).*abovep;
cc(15,ix) = cc(15,ix) + adt*(relabel.dM12deta2 - relabel.dM21deta2).*abovep;
cc(16,ix) = cc(16,ix) + adt*(relabel.dM12dq1   - relabel.dM21dq1  ).*abovep;
cc(17,ix) = cc(17,ix) + adt*(relabel.dM12dq2   - relabel.dM21dq2  ).*abovep;


% sigma equation

ix = 9:9:9*nz;
cc(4 ,ix) = cc(4 ,ix) - m1.*eos.drdetap1.*belowp;
cc(5 ,ix) = cc(5 ,ix) - m2.*eos.drdetap2.*belowp;
cc(6 ,ix) = cc(6 ,ix) - m1.*eos.drdqp1  .*belowp;
cc(7 ,ix) = cc(7 ,ix) - m2.*eos.drdqp2  .*belowp;
cc(8 ,ix) = cc(8 ,ix) + 1./eos.rho1;
cc(9 ,ix) = cc(9 ,ix) + 1./eos.rho2;
cc(10,ix) = - m1.*eos.drdp1 - m2.*eos.drdp2; 
cc(13,ix) = cc(13,ix) - m1.*eos.drdetap1.*abovep;
cc(14,ix) = cc(14,ix) - m2.*eos.drdetap2.*abovep;
cc(15,ix) = cc(15,ix) - m1.*eos.drdqp1  .*abovep;
cc(16,ix) = cc(16,ix) - m2.*eos.drdqp2  .*abovep;

