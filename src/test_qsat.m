% Test the calculation of qsat

p = 97000;

T = 297.0;

%for eta = 1110:0.1:1112
for eta = 1111.001
    
    eta
    [qsat,dqsatdeta ] = find_qsatl(p,eta,T,constants.therm)
    
end