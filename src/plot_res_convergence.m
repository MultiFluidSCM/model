% Plot of residuals to look at convergence
    
figure(30)
fz = 16;
irange = [1:qn_iter];

subplot(3,3,1)
norm1 = maxres1w(1);
norm2 = maxres2w(1);
plot(irange,maxres1w(irange)/norm1,'bo',irange,maxres2w(irange)/norm2,'ro')
title('Res w')
xlim([0,qn_iter_max])
set(gca,'fontsize',fz)

subplot(3,3,2)
norm1 = maxres1m(1);
norm2 = maxres2m(1);
plot(irange,maxres1m(irange)/norm1,'bo',irange,maxres2m(irange)/norm2,'ro')
title('Res m')
xlim([0,qn_iter_max])
set(gca,'fontsize',fz)

subplot(3,3,3)
norm1 = maxres1eta(1);
norm2 = maxres2eta(1);
plot(irange,maxres1eta(irange)/norm1,'bo',irange,maxres2eta(irange)/norm2,'ro')
title('Res eta')
xlim([0,qn_iter_max])
set(gca,'fontsize',fz)

subplot(3,3,4)
norm1 = maxres1q(1);
norm2 = maxres2q(1);
plot(irange,maxres1q(irange)/norm1,'bo',irange,maxres2q(irange)/norm2,'ro')
title('Res q')
xlim([0,qn_iter_max])
set(gca,'fontsize',fz)

subplot(3,3,5)
norm1 = maxres1u(1);
norm2 = maxres2u(1);
norm3 = maxres1v(1);
norm4 = maxres2v(1);
plot(irange,maxres1u(irange)/norm1,'bo',irange,maxres2u(irange)/norm2,'ro',...
     irange,maxres1v(irange)/norm3,'b+',irange,maxres2v(irange)/norm4,'r+')
title('Res u, v')
xlim([0,qn_iter_max])
set(gca,'fontsize',fz)

subplot(3,3,6)
norm1 = maxres1eosetap(1);
norm2 = maxres2eosetap(1);
norm3 = maxres1eoseta(1);
norm4 = maxres2eoseta(1);
norm5 = maxressigma(1);
plot(irange,maxres1eosetap(irange)/norm1,'bo',irange,maxres2eosetap(irange)/norm2,'ro',...
     irange,maxres1eoseta(irange)/norm3,'b+',irange,maxres2eoseta(irange)/norm4,'r+',...
     irange,maxressigma(irange)/norm5,'ko')
title('Res eos')
xlim([0,qn_iter_max])
set(gca,'fontsize',fz)

subplot(3,3,7)
norm1 = maxres1tke(1);
norm2 = maxres2tke(1);
plot(irange,maxres1tke(irange)/norm1,'bo',irange,maxres2tke(irange)/norm2,'ro')
title('Res tke')
xlim([0,qn_iter_max])
set(gca,'fontsize',fz)

subplot(3,3,8)
norm1 = maxres1Veta(1);
norm2 = maxres2Veta(1);
plot(irange,maxres1Veta(irange)/norm1,'bo',irange,maxres2Veta(irange)/norm2,'ro')
title('Res Veta')
xlim([0,qn_iter_max])
set(gca,'fontsize',fz)

subplot(3,3,9)
norm1 = maxres1Vq(1);
norm2 = maxres2Vq(1);
plot(irange,maxres1Vq(irange)/norm1,'bo',irange,maxres2Vq(irange)/norm2,'ro')
title('Res Vq')
xlim([0,qn_iter_max])
set(gca,'fontsize',fz)
