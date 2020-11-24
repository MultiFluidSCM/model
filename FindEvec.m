function [ ] = FindEvec( cc, grid )
%Ndiagsolvex Solve an N-diagonal linear system
%   This version provides a sanity check on the version Ndiagsolveb
%   by converting the compact representation of the matrix to a full
%   explicit representation and using Matlabs inbuilt linear solver
% The coefficients in row i are cc(1:2*p+1,i)

[Ndiag,N] = size(cc);
p = (Ndiag-1)/2;
d = p + 1;

% Construct explicit matrix
A = zeros(N,N);
for i = 1:N
    j1 = max(d+1-i,1);
    j2 = min(d+N-i,Ndiag);
    k1 = max(i-p,1);
    k2 = min(i+p,N);
    A(i,k1:k2) = cc(j1:j2,i);
end
disp('Condition number of matrix:')
cond(A)

[V,D] = eig(A);
l0 = diag(D);
[lambda,idx] = sort(abs(l0));

figure(16)
plot(lambda)
title('Eigenvalue spectrum')
lambda(1:10)
pause

nz = grid.nz;

ix = 8;
disp(['Evals ',num2str(ix),' and ',num2str(N)])
l0(idx(ix))
l0(idx(N))

xx = V(:,idx(ix));
inc_w1   = xx(1:9:9*nz+1);
inc_w2   = xx(2:9:9*nz+2);
inc_eta1 = xx(3:9:9*nz+3);
inc_eta2 = xx(4:9:9*nz+4);
inc_q1   = xx(5:9:9*nz+5);
inc_q2   = xx(6:9:9*nz+6);
inc_m1   = xx(7:9:9*nz-2);
inc_m2   = xx(8:9:9*nz-1);
inc_p    = xx(9:9:9*nz);

figure(17)
prange = 1:nz;
fs = 15;
subplot(2,3,1)
plot(real(inc_p(prange)),grid.zp(prange),'b',imag(inc_p(prange)),grid.zp(prange),'b:')
title(['p evec mode ',num2str(ix)])
set(gca,'fontsize',fs)
subplot(2,3,2)
plot(real(inc_w1(prange)),grid.zw(prange),'b',imag(inc_w1(prange)),grid.zw(prange),'b:',...
     real(inc_w2(prange)),grid.zw(prange),'r',imag(inc_w2(prange)),grid.zw(prange),'r:')
title('w evec')
set(gca,'fontsize',fs)
subplot(2,3,3)
plot(real(inc_m1(prange)),grid.zp(prange),'b',imag(inc_m1(prange)),grid.zp(prange),'b:',...
     real(inc_m2(prange)),grid.zp(prange),'r',imag(inc_m2(prange)),grid.zp(prange),'r:')
title('m evec')
set(gca,'fontsize',fs)
subplot(2,3,4)
plot(real(inc_eta1(prange)),grid.zw(prange),'b',imag(inc_eta1(prange)),grid.zw(prange),'b:',...
     real(inc_eta2(prange)),grid.zw(prange),'r',imag(inc_eta2(prange)),grid.zw(prange),'r:')
title('eta evec')
set(gca,'fontsize',fs)
subplot(2,3,5)
plot(real(inc_q1(prange)),grid.zw(prange),'b',imag(inc_q1(prange)),grid.zw(prange),'b:',...
     real(inc_q2(prange)),grid.zw(prange),'r',imag(inc_q2(prange)),grid.zw(prange),'r:')
title('q evec')
set(gca,'fontsize',fs)
pause
figure(1)
pause



end

