% Sample fluid 2 profile and plot buoyancy of various parcels
% lifted adiabatically. Assume environment p and fluid 1 rho
% remain unchanged as parcel is lifted.

figure(19)

% First plot actual buoyancy profile
plot(work.buoy,zunitsw,'k','linewidth',1.5)
set(gca,'FontSize',fs)
ylim([0,zplottop])
title('Updraft buoyancy ')
xlabel('b')
ylabel(labelz)
hold on


% List of starting levels (w-levels)
kstart = [2,30,35];
% kstart = [35,35,35,35]; % Use this line with the addition of pertubations below
% kstart = [43];
count = 0;

% Loop over starting levels
for k0 = kstart
    
    count = count + 1;
    
    % Parcel eta and q will be conserved
    eta_parcel = state_new.fluid(2).eta(k0);
    q_parcel   = state_new.fluid(2).q(k0);
%     if count == 2 | count == 4
%         eta_parcel = eta_parcel + sqrt(state_new.fluid(2).vareta(k0));
%     end
%     if count > 2
%         q_parcel   = q_parcel   + sqrt(state_new.fluid(2).varq(k0));
%     end
    
    % Lift parcel, calculate density, hence buoyancy
    % nk = 30;
    nk = 70 - k0;
    % nk = 1;
    buoy_parcel = zeros(1,nk);
    for ik = 1:nk
        
        k = k0 + (ik - 1);
        
        % Compute density of parcel at level k
        % First we need to solve for T
        p_parcel = eos.pbar(k);
        T_parcel = state_new.fluid(2).Tw(k);
        for iter = 1:5
            [g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(p_parcel, ...
                                                           T_parcel, ...
                                                           q_parcel,  ...
                                                           constants.therm);
            res = gt + eta_parcel
            T_inc = -res/gtt;
            T_parcel = T_parcel + T_inc;
        end
        rho_parcel = 1/gp;
        buoy_parcel(ik) = constants.phys.gravity*sigma1(k)*(eos.rhow1(k) - rho_parcel)/rho_parcel;
        is_liquid(ik) = 1 - a < q_parcel - 1e-7;;
        disp('check T')
        [state_new.fluid(2).Tw(k0),T_parcel]
        disp('check rho')
        [eos.rhow2(k0),rho_parcel]
        disp('liquid')
        (q_parcel + a - 1)/a
        %pause
        
    end
    
    % Plot buoyancy of lifted parcel
    plot(buoy_parcel,zunitsw(k0:k0+nk-1),'b-o')
    for ik = 1:nk
        if is_liquid(ik)
            plot([buoy_parcel(ik)],[zunitsw(k0+ik-1)],'r-o')
        end
    end
    
end

hold off

