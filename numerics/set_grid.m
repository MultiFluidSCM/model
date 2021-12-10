function grid = set_grid(grid_settings)
% Set up computational domain and grid

% Number of levels
nz = grid_settings.nz;
nzp = nz + 1;

% Domain range
zbot = grid_settings.zbot;
ztop = grid_settings.ztop;

% Stretch factor for grid
% Ratio of top layer thickness to lowest layer thickness
totstretch = grid_settings.totstretch;
% Ratio of adjacent layer thicknesses
stretch = totstretch^(1/(nz-1));
rs = sqrt(stretch);

% Set up grid
if stretch == 1.0
    hdz0 = ztop/(2*nz);
else
    hdz0 = ztop*(1.0 - rs)/(1.0 - rs^(2*nz));
end
hdz = hdz0;

grid.zw(1) = zbot;
for k = 1:nz
    grid.zp(k) = grid.zw(k) + hdz;
    hdz = hdz*rs;
    grid.zw(k+1) = grid.zp(k) + hdz;
    hdz = hdz*rs;
end

% Grid in km for plotting
grid.zwkm = grid.zw*1e-3;
grid.zpkm = grid.zp*1e-3;

% Grid intervals
grid.dzp = grid.zw(2:nzp) - grid.zw(1:nz);
grid.dzw(1) = grid.zp(1) - grid.zw(1);
grid.dzw(2:nz) = grid.zp(2:nz) - grid.zp(1:nz-1);
grid.dzw(nzp) = grid.zw(nzp) - grid.zp(nz);
disp(['Min dzp ',num2str(min(grid.dzp)),' Max dzp ',num2str(max(grid.dzp))])
%pause

% Interpolation coefficients
% From p to w levels
grid.abovew(1) = 1.0;
grid.beloww(1) = 0.0;
for k = 2:nz
    grid.abovew(k) = (grid.zw(k) - grid.zp(k-1))/(grid.zp(k) - grid.zp(k-1));
    grid.beloww(k) = 1.0 - grid.abovew(k);
end
grid.abovew(nzp) = 0.0;
grid.beloww(nzp) = 1.0;
% and for extrapolation to the bottom and top boundaries
grid.extrapb1 = (grid.zp(2) - grid.zw(1))/grid.dzw(2);
grid.extrapb2 = 1.0 - grid.extrapb1;
grid.extraptnz = (grid.zw(nzp) - grid.zp(nz-1))/grid.dzw(nz);
grid.extraptnzm = 1.0 - grid.extraptnz;
% From w to p levels
for k = 1:nz
    grid.abovep(k) = (grid.zp(k) - grid.zw(k))/(grid.zw(k+1) - grid.zw(k));
    grid.belowp(k) = 1.0 - grid.abovep(k);
end
% Remapping coefficients
% From p to w levels
grid.abover = grid.beloww;
grid.belowr = grid.abovew;
grid.abover(1) = 1.0;
grid.belowr(1) = 0.0;
grid.abover(nzp) = 0.0;
grid.belowr(nzp) = 1.0;
% From w to p levels
grid.aboves = grid.belowp;
grid.belows = grid.abovep;

% Remember inputs for later
grid.nz = nz;
grid.nzp = nzp;
grid.zbot = zbot;
grid.ztop = ztop;
grid.totstretch = totstretch;

% test_product_rule(grid)

end