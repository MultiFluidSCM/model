% One-off test of gibbs and derivatives

% p = 9.6998297e4;
% T = 298.636932;
% q = 0.0151957;
% a = 1 - q;
% 
% dp = 0.0377;
% dT = -0.369393;
% da = -0.000702;
% 
% 
% [g,gp,gt,ga,gpp,gpt,gpa,gtt,gta,gaa,gtwv] = ...
%     gibbsav(p,T,a,constants.therm);
% RH = relhum(p,T,1-a,constants.therm)
% 
% [gx,gpx,gtx,gax,gppx,gptx,gpax,gttx,gtax,gaax,gtwvx] = ...
%     gibbsav(p+dp,T,a+da,constants.therm);
% RHx = relhum(p+dp,T,1-a-da,constants.therm)
% 
% 
% disp(' deta ')
% - (gtx - gt)
% - gta*da - gpt*dp - gtt*dT
% disp(' dalpha ')
% gpx - gp
% gpa*da + gpp*dp + gpt*dT
% disp('drho')
% 1/gpx - 1/gp
% 
% gw = g - a*ga
% etaw = -(gt - a*gta)
% hw = gw + T*etaw
% 
% %gtwv
% 
% a = 1e-6;
% [g,gp,gt,ga,gpp,gpt,gpa,gtt,gta,gaa,gtwv] = ...
%     gibbsav(p,T,a,constants.therm);
% 
% gw_hat = g - a*ga
% etaw_hat = -(gt - a*gta)
% hw_hat = gw_hat + T*etaw_hat
% 
% %gtwv
% 
% eta_per_dm = (gw_hat - gw + T*etaw_hat)/T
% 
% eta_per_dm = (hw_hat - hw + T*etaw)/T
% 
% disp('check')
% gw     + T*etaw     - hw
% gw_hat + T*etaw_hat - hw_hat

disp(' ')
disp(' ')
% Starting values
p = 97000;
T = 296.4183
q = 0.01497;
a = 1 - q;
rho = 1.129954604805245;
V = 0.01*22.07;
M = rho*V;
Mv = q*M;
[g,gp,gt,ga,gpp,gpt,gpa,gtt,gta,gaa,gtwv] = ...
    gibbsav(p,T,a,constants.therm);
RH = relhum(p,T,1-a,constants.therm);
disp('rho (bulk)')
[rho,1/gp]
disp('eta (bulk)')
-gt



eta = -gt;
etaw = -(gt - a*gta)
alphaw = (gp - a*gpa);

% Initial energy (bulk)
e = (g - T*gt - p*gp)
E = e*M
% and entropy
N = -gt*M

% Specific energy of water
ew = (g - T*gt - p*gp) - a*(ga -T*gta - p*gpa)
% constants.therm.Cvv*T + constants.therm.L00
% and dry air
ed = constants.therm.Cvd*T

% Added mass and volume of water
dM = 0.00012;
dV = alphaw*dM;

% Added energy
dE = ew*dM
Etot = E + dE
% Added entropy
dN = etaw*dM

% Net quantities in initial state
rho_net = (M + dM)/(V + dV)
e_net = Etot/(M + dM)
n_net = (N + dN)/(M + dM)

% Work done
W = p*dV


% Now compress into volume V
% New state
M_new = M + dM;
rho_new = M_new/V
Mv_new = Mv + dM;
q_new = Mv_new/M_new;
a_new = 1 - q_new;
e_new = (Etot + W)/(M + dM)
n_new = (N + dN)/(M + dM)  % Naive calculation?

disp(' ')
disp(' using energy ')
p_new = p;
T_new = T;
for iter = 1:5
    [g,gp,gt,ga,gpp,gpt,gpa,gtt,gta,gaa,gtwv] = ...
        gibbsav(p_new,T_new,a,constants.therm);
    r1 = gp - 1/rho_new;
    r2 = (g - T_new*gt - p_new*gp) - e_new;
    A11 = gpp;
    A12 = gpt;
    A21 = (- T_new*gpt - p_new*gpp);
    A22 = (- T_new*gtt - p_new*gpt);
    rdet = 1/(A11*A22 - A12*A21);
    dp = rdet*(-A22*r1 + A12*r2);
    dT = rdet*( A21*r1 - A11*r2);
    p_new = p_new + dp;
    T_new = T_new + dT;
end
disp('new rho')
[rho_new,1/gp]
disp('new eta')
[n_new, -gt]
disp('New p')
p_new
disp('New T')
T_new
e_new = (g - T_new*gt - p_new*gp)
E_new = e_new*(M + dM)


disp(' ')
disp(' using entropy ')
eta_new = n_new
p_new = p;
T_new = T;
for iter = 1:5
    [g,gp,gt,ga,gpp,gpt,gpa,gtt,gta,gaa,gtwv] = ...
        gibbsav(p_new,T_new,a,constants.therm);
    r1 = gp - 1/rho_new;
    r2 = gt + eta_new;
    rdet = 1/(gpp*gtt - gpt*gpt);
    dp = rdet*(-gtt*r1 + gpt*r2);
    dT = rdet*( gpt*r1 - gpp*r2);
    p_new = p_new + dp;
    T_new = T_new + dT;
end
disp('new rho')
[rho_new,1/gp]
disp('new eta')
[eta_new, -gt]
disp('New p')
p_new
disp('New T')
T_new
e_new = (g - T_new*gt - p_new*gp)
E_new = e_new*(M + dM)


% % Now return adiabtatically to original pressure
% 
% p_fin = p;
% T_fin = T;
% for iter = 1:5
%     
%     [g,gp,gt,ga,gpp,gpt,gpa,gtt,gta,gaa,gtwv] = ...
%         gibbsav(p_fin,T_fin,a,constants.therm);
% 
%     r2 = gt + eta_new;
%     dT = -r2/gtt;
%     T_fin = T_fin + dT;
% 
% end
% 
% disp('fin rho')
% 1/gp
% disp('fin eta')
% -gt
% disp('fin p')
% p_fin
% disp('fin T')
% T_fin

% Now return adiabtatically to original volume

disp(' ')
p_fin = p_new;
T_fin = T_new;
eta_fin = eta_new;
rho_fin = (M + dM)/(V + dV);
for iter = 1:5
    
    [g,gp,gt,ga,gpp,gpt,gpa,gtt,gta,gaa,gtwv] = ...
        gibbsav(p_fin,T_fin,a,constants.therm);

    r1 = gp - 1/rho_fin;
    r2 = gt + eta_fin;
    rdet = 1/(gpp*gtt - gpt*gpt);
    dp = rdet*(-gtt*r1 + gpt*r2);
    dT = rdet*( gpt*r1 - gpp*r2);
    p_fin = p_fin + dp;
    T_fin = T_fin + dT;

end

disp('fin rho')
1/gp
disp('fin eta')
-gt
disp('fin p')
p_fin
disp('fin T')
T_fin

e_fin = (g - T_fin*gt - p_fin*gp)
E_fin = e_fin*(M + dM)



