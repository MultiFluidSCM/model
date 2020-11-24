function [  ] = test_product_rule( grid )
% Test discrete product rule for derivatives on w levels

nz = numel(grid.zp);
nzp = nz+1;

grid.abovew + grid.abover
pause
grid.beloww + grid.belowr
pause

% Set up functions to work with
f1 = grid.zp.*grid.zp;
f2 = exp( 0.00001*grid.zp);
f3 = exp( 0.00001*grid.zw);
f1bar = weight_to_w(grid,f1);

% Product rule for d/dz (f1*f2)

dfull(1) = 0;
dfull(2:nz) = f1(2:nz).*f2(2:nz) - f1(1:nz-1).*f2(1:nz-1);
dfull(nzp) = 0

term1(1) = 0;
term1(2:nz) = f1bar(2:nz).*(f2(2:nz) - f2(1:nz-1));
term1(nzp) = 0;

% f2bar(1) = 0;
% f2bar = grid.abovew(2:nz).*f2(2:nz) + grid.beloww(2:nz).*f2(1:nz-1);
% f2bar(nzp) = 0;

f2bar(1:nz) = grid.abovew(1:nz).*f2;
f2bar(nzp) = 0;
f2bar(2:nzp) = f2bar(2:nzp) + grid.beloww(2:nzp).*f2;

term2(1) = 0;
term2(2:nz) = f2bar(2:nz).*(f1(2:nz) - f1(1:nz-1));
term2(nzp) = 0;

term1 + term2
pause

% Product rule for d/dz( f1bar *f3 )

dfull = (f1bar(2:nzp).*f3(2:nzp) - f1bar(1:nz).*f3(1:nz))

term1 = f1.*(f3(2:nzp) - f3(1:nz));

df1(1) = 0;
df1(2:nz) = (f1(2:nz) - f1(1:nz-1));
df1(nzp) = 0;
temp = f3.*df1;
term2 = grid.beloww(2:nzp).*temp(2:nzp) + grid.abovew(1:nz).*temp(1:nz);

term1 + term2
pause



end

