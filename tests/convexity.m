% Test convexity of Gibbs function

format long

constants.therm = set_therm_const();

p = 1e5;
T = 300;
dT = 0.01;
dp = 1;
q = 0.01;
dq = 0.001;

[g,gp,gt,gw,gpp,gpt,gtt,gpw,gtw,gww,a] = gibbs(p   ,T   ,q   ,constants.therm);
g_air = g - q*gw;
g_water = gw + g_air;
c1 = constants.therm.epsilon*(1 - q);
c2 = q;
p_air = p*c1/(c1 + c2);
p_water = p*c2/(c1 + c2);


[G,Gp,Gt,Gw,Gpp,Gpt,Gtt,Gpw,Gtw,Gww,a] = gibbs(p   ,T+dT,q   ,constants.therm);
G_air = G - q*Gw;
G_water = Gw + G_air;
Gt_air = Gt - q*Gtw;
Gt_water = Gtw + Gt_air;
alpha_air = Gp/(1 - q);
alpha_water = Gp/q;
c1 = constants.therm.epsilon*(1 - q);
c2 = q;
P_air = p*c1/(c1 + c2);
P_water = p*c2/(c1 + c2);

curvT = G - g - Gt*dT
curvT_air   = G_air   - g_air   - alpha_air  *(P_air   - p_air)   - Gt_air*dT
curvT_water = G_water - g_water - alpha_water*(P_water - p_water) - Gt_water*dT

% ---

[G,Gp,Gt,Gw,Gpp,Gpt,Gtt,Gpw,Gtw,Gww,a] = gibbs(p+dp,T   ,q   ,constants.therm);
G_air = G - q*Gw;
G_water = Gw + G_air;
Gt_air = Gt - q*Gtw;
Gt_water = Gtw + Gt_air;
alpha_air = Gp/(1 - q);
alpha_water = Gp/q;
c1 = constants.therm.epsilon*(1 - q);
c2 = q;
P_air = (p + dp)*c1/(c1 + c2);
P_water = (p + dp)*c2/(c1 + c2);

curvp = G - g - Gp*dp
curvp_air   = G_air   - g_air   - alpha_air  *(P_air   - p_air)   - Gt_air*dT
curvp_water = G_water - g_water - alpha_water*(P_water - p_water) - Gt_water*dT

% ---

[G,Gp,Gt,Gw,Gpp,Gpt,Gtt,Gpw,Gtw,Gww,a] = gibbs(p+dp,T+dT,q   ,constants.therm);
G_air = G - q*Gw;
G_water = Gw + G_air;
Gt_air = Gt - q*Gtw;
Gt_water = Gtw + Gt_air;
alpha_air = Gp/(1 - q);
alpha_water = Gp/q;
c1 = constants.therm.epsilon*(1 - q);
c2 = q;
P_air = (p + dp)*c1/(c1 + c2);
P_water = (p + dp)*c2/(c1 + c2);

curvpT = G - g - Gp*dp -Gt*dT
curvpT_air   = G_air   - g_air   - alpha_air  *(P_air   - p_air)   - Gt_air*dT
curvpT_water = G_water - g_water - alpha_water*(P_water - p_water) - Gt_water*dT


format short
