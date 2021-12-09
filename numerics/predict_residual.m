% Predict residual in u1 equation

disp(' ')
disp('Predicted residual')

% Tendency term
dd = zeros(3,nz);
dd(2,ix) = dd(2,ix) + m1;
yy = state_err.fluid(1).u;
rr1 = - Ndiagmult(dd,yy);
disp(['Tendency    ' num2str(rr1(krange))])

% Transport terms
dd = zeros(3,nz);
dd(1,ix) = dd(1,ix) + adt*m1.*work.dAu1dub;   
dd(2,ix) = dd(2,ix) + adt*m1.*work.dAu1duc;   
dd(3,ix) = dd(3,ix) + adt*m1.*work.dAu1dua;
yy = state_err.fluid(1).u;
rr2 = - Ndiagmult(dd,yy);
disp(['Transport   ' num2str(rr2(krange))])

% Diffusion terms
dd = zeros(3,nz);
dd(1,ix) = dd(1,ix) - adt*work.dDu1dub(1:nz)./dzp;
dd(2,ix) = dd(2,ix) - adt*(work.dDu1dua(1:nz) - work.dDu1dub(2:nzp))./dzp;
dd(3,ix) = dd(3,ix) + adt*work.dDu1dua(2:nzp)./dzp;
yy = state_err.fluid(1).u;
rr3 = - Ndiagmult(dd,yy);
disp(['Diffusion   ' num2str(rr3(krange))])

% Relabelling terms
rr4 = - (-adt*M12.*(work.duhat12du1 - 1) + adt*M21.*(work.duhat21du1 - 1)).*state_err.fluid(1).u ...
      - (-adt*M12.*(work.duhat12du2)     + adt*M21.*(work.duhat21du2)    ).*state_err.fluid(2).u;
disp(['E minus D   ' num2str(rr4(krange))])


rr = rr1 + rr2 + rr3 + rr4;
disp(['Total       ',num2str(rr(krange))])

disp(' ')
disp('Actual residual')
rr1 = -m1.*state_err.fluid(1).u;
disp(['Tendency    ' num2str(rr1(krange))])
rr2 = (tend.fluid(1).u.transport - tend_tru.fluid(1).u.transport)*adt.*m1;
disp(['Transport   ' num2str(rr2(krange))])
rr3 = (tend.fluid(1).u.diffuse   - tend_tru.fluid(1).u.diffuse  )*adt.*m1;
disp(['Diffusion   ' num2str(rr3(krange))])
rr4 = (tend.fluid(1).u.relabel   - tend_tru.fluid(1).u.relabel  )*adt.*m1;
disp(['E minus D   ' num2str(rr4(krange))])
rr = rr1 + rr2 + rr3 + rr4;
disp(['Total       ',num2str(rr(krange))])

pause