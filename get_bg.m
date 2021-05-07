function [alpha_b,beta_b,gamma_b,eta] = get_bg(x,T,Tc,rho,rhoc,p,pc,M)
% calculate the background coefs

%% preperation

% convert units
M    = M*1000;      % from kg/mol to g/mol
rho  = rho/M;       % from kg/m3 to mol/L
rhoc = rhoc/M;      % from kg/m3 to mol/L
pc   = pc*1e-6;     % from Pa to MPa
p    = p*1e-6;      % from Pa to MPa

% constants
R  = 8314.51;
M1 = 44.01;
M2 = 30.069;
Tc1 = 304.136;
Tc2 = 305.322;
Pc1 = 7.3773;
Pc2 = 4.8718;
rhoc1 = 10.625;
rhoc2 = 6.8701;

% pure conponent properties, use original units 
[eta_1,lamda_01,lamda_ex1] = pure_conp(T,rho,1);
[eta_2,lamda_02,lamda_ex2] = pure_conp(T,rho,2);

%% calculate
% the compressibility factor
Zcx = pc/R/rhoc/Tc*1e6;
% Eq.(33) in Kiselev, Huber. Fluid phase equilibria, 1998, 142.
Gammacx  = Tc^(1/6)*Zcx^5*sqrt(M)/pc^(2/3);
% Eq.(39)
D0 = 1.01325e-8*T^1.75*sqrt((M1+M2)/M1/M2)/p/( 26.90^(1/3) + 44.88^(1/3) )^2;
% Eq.(31)
alpha_coef = [0,0,0,0,0,1.19792156449951e-05,1.12168634915470e-05,0,0,0,0,-4.53763235100211e-06,-3.29427570391861e-06,0,0,0,0,5.64431267032387e-07,1.16814090337172e-07];
SUM = 0;
for ind = [1,2,3,4,5,6]
    SUM = SUM + ( alpha_coef(3*ind) + alpha_coef(3*ind+1)*x )*(rho/rhoc)^(ind+1);
end
alpha_ex = x*(1-x)*R^(-7/6)/Tc/Gammacx*SUM;
% Eq.(32)
beta_coef = [-1.93675085936116,0,-3.29858643973719e-06,0,0,2.50317378116377e-06,0,0,0,0,0,-2.91075290773087e-07,0,0,0,0,0,1.32106664277325e-08,0];
SUM = 0;
for ind = [1,2,3,4,5,6]
    SUM = SUM + ( beta_coef(3*ind) + beta_coef(3*ind+1)*x )*(rho/rhoc)^(ind+1);
end
beta_ex = x*(1-x)*R^(-1/6)/Tc/Gammacx*SUM;
% Eq.(36); the current dimensionSets are correct
alpha_0 = rho*D0*M^2/R/T*x*(1-x)*( 0.504374247884148 - 0.334397719244566*(rho/rhoc)^3 + 0.123793083165651*(rho/rhoc)^5 - 0.0130863246445822*(rho/rhoc)^7 );
% Eq.(37); must divide M to consist the dimensionSet
beta_0  = R/M*alpha_0*(beta_coef(1)+beta_coef(2)*x-log(x/(1-x)));
% Eq.(29) and (30)
alpha_b = alpha_0 + alpha_ex;
beta_b  = beta_0  + beta_ex;

% Eq.(43)
S1 = 1.5*194.7;
S2 = 1.5*184.5;
S12 = sqrt(S1*S2);
% Eq.(41)
A12 = (1+sqrt(eta_1/eta_2*(M2/M1)^(3/4)*(T+S1)/(T+S2)))^2*(T+S12)/(T+S1)/4;
A21 = (1+sqrt(eta_2/eta_1*(M1/M2)^(3/4)*(T+S2)/(T+S1)))^2*(T+S12)/(T+S2)/4;
% Eq.(40)
lamda_0 = (1-x)*lamda_01/(1-x+x*A12) + x*lamda_02/(x+(1-x)*A21);
% Eq.(35)
Zc1 = Pc1/R/rhoc1/Tc1*1e6;
Zc2 = Pc2/R/rhoc2/Tc2*1e6;
lamda_ex = lamda_ex1*(1-x)*(Tc1^(1/6)*Zc1^5*sqrt(M1)/Gammacx/Pc1^(2/3)) + lamda_ex2*x*(Tc2^(1/6)*Zc2^5*sqrt(M2)/Gammacx/Pc2^(2/3));
% Eq.(34)
gamma_coef = [-1.17562e-6,-5.14420e-7,5.54845e-6,0,-6.72675e-6,2.60059e-6,0,0,1.74581e-6,-1.24763e-6,0,0];
SUM = 0;
for ind = [1,2,3,4,5,6]
    SUM = SUM + ( gamma_coef(2*ind-1) + gamma_coef(2*ind)*x )*(rho/rhoc)^ind;
end
gamma_b = lamda_0 + lamda_ex + T*beta_b^2/alpha_b + x*(1-x)*R^(5/6)/Gammacx*SUM;

% Eq.(85)
eta = (eta_1*(1-x)*Tc1^(1/6)/sqrt(M1)/Pc1^(2/3) + eta_2*x*Tc2^(1/6)/sqrt(M2)/Pc2^(2/3))*sqrt(M)*pc^(2/3)/Tc^(1/6);

end





    



