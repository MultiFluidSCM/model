function zstar = findzstar(zw,w2,Tw,theta_rho1,theta_rho2,option_zstar)

% Estimate the inversion height zstar

% Threshold for options 1 and 4
thresh = 0.02*max(w2);
% Domain top
ztop = max(zw);

% Option 1
% Height at which w2 becomes zero (or < some threshold)
k = 2; 
while (w2(k) > thresh)
    k = k + 1;
end
zstar1 = zw(k);
if (k > 2)
    km = k - 1;
    zstar1 = zw(km) - (w2(km) - thresh)*(zw(k)-zw(km))/(w2(k) - w2(km));
end


if option_zstar == 2

    % Option 2
    % Lowest height at which extrapolated w2 becomes zero
    % Linear extrapolation
    k = 2;
    zstar2 = min(2*zstar1,ztop);
    while (w2(k) > 0 & zw(k) < zstar2)
        k = k + 1;
        xdw = w2(k) - w2(k-1);
        xdz = zw(k) - zw(k-1);
        if (xdw < 0)
            z0 = zw(k) - w2(k)*xdz/xdw;
            zstar2 = min(zstar2,z0);
        end
    end

    
elseif option_zstar == 3

    % Option 3
    % Lowest height at which extrapolated w2 becomes zero
    % Fit to (z* - z)^0.5
    k = 2;
    zstar3 = min(2*zstar1,ztop);
    while (w2(k) > 0 & zw(k) < zstar3)
        k = k + 1;
        xdw = w2(k) - w2(k-1);
        xdz = zw(k) - zw(k-1);
        if (xdw < 0)
            wsqa = w2(k)*w2(k);
            wsqb = w2(k-1)*w2(k-1);
            z0 = (wsqa*zw(k-1) - wsqb*zw(k))/(wsqa - wsqb);
            if z0 < zstar3 & z0 > zw(k)
                zstar3 = z0;
            end
        end
    end
    
    
elseif option_zstar == 4

    % Look for maximum dT/dz
    % First find which interval it is in
%     k = 2;
    kmxdtdz = 2;
    mxdtdz = -0.01;
    for k=2:length(zw)-1
        dtdz = (Tw(k+1) - Tw(k))/(zw(k+1)-zw(k));
        if (dtdz > mxdtdz)
            kmxdtdz = k;
            mxdtdz = dtdz;
        end
    end
%     while (w2(k) > thresh)
%         k = k + 1;
%         dtdz = (theta1(k+1) - theta1(k))/(zw(k+1)-zw(k));
%         if (dtdz > mxdtdz)
%             kmxdtdz = k;
%             mxdtdz = dtdz;
%         end
%     end
    % Now find the exact location by a cubic fit
    % Lagrange cubic
%     zzm = zw(kmxdtdz-1);
%     zz0 = zw(kmxdtdz  );
%     zzp = zw(kmxdtdz+1);
%     zzpp = zw(kmxdtdz+2);
%     fac1 = Tw(kmxdtdz-1)/((zzm-zz0)*(zzm-zzp)*(zzm-zzpp));
%     fac2 = Tw(kmxdtdz  )/((zz0-zzm)*(zz0-zzp)*(zz0-zzpp));
%     fac3 = Tw(kmxdtdz+1)/((zzp-zzm)*(zzp-zz0)*(zzp-zzpp));
%     fac4 = Tw(kmxdtdz+2)/((zzpp-zzm)*(zzpp-zz0)*(zzpp-zzp));
%     zstar4 = (fac1*(zz0+zzp+zzpp) + fac2*(zzm+zzp+zzpp) ...
%            + fac3*(zzm+zz0+zzpp) + fac4*(zzm+zz0+zzp )) ...
%           / (3*(fac1 + fac2 + fac3 + fac4));
%     zstar4 = max([zstar4,zz0]);
%     zstar4 = min([zstar4,zzp]);
    
    % Hermite cubic
    zzm = zw(kmxdtdz-1);
    zz0 = zw(kmxdtdz  );
    zzp = zw(kmxdtdz+1);
    zzpp = zw(kmxdtdz+2);
    T0 = Tw(kmxdtdz  );
    Tp = Tw(kmxdtdz+1);
    dz = zzp - zz0;
    dTp = (Tw(kmxdtdz+2) - Tw(kmxdtdz  ))/(zzpp - zz0 );
    dT0 = (Tw(kmxdtdz+1) - Tw(kmxdtdz-1))/(zzp  - zzm );
    zhat = -(3*(Tp - T0) - dz*(dTp + 2*dT0))/(3*(dz*(dTp + dT0) -2*(Tp - T0)));
    zstar4 = zz0 + dz*zhat;
    % Plot to check
    % plot_cubic_fit
    
    % This option has issues
    % 1. It can jump to a level remote from the true BL top (e.g. ARM case)
    % 2. It tends to pick out the middle of the layer so zstar ascends in
    % steps

elseif option_zstar == 5

    % Where does theta equal surface theta?
    % Search for first theta greater than thetas
    % Work with density potential temperature for consistency
    k=2;
    theta_rho_s = theta_rho2(1);
    while theta_rho1(k) <= theta_rho_s
        k = k + 1;
    end
    
    zstar5 = zw(k-1) + (zw(k) - zw(k-1))*(theta_rho_s - theta_rho1(k-1))/(theta_rho1(k) - theta_rho1(k-1));
    zstar5 = max(zstar5,50.0);
    
end



if option_zstar == 1
    zstar = zstar1;
elseif option_zstar == 2
    zstar = zstar2;
elseif option_zstar == 3
    zstar = zstar3;
elseif option_zstar == 4
    zstar = zstar4;
else
    zstar = zstar5;       
end


end

    
    
    