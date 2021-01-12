[x,y] = meshgrid(0:0.01:2, 0:0.01:3);

z = 4*log(x)+ log(y) - 6*x - y;

contour(x,y,z)