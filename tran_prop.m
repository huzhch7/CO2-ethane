function [eta,lamda,D,kT] = tran_prop(x,T,Tc,rho,rhoc,p,pc,cp,mu_x,mu_T,ALPHA,KAPPA)

M1 = 0.04401;
M2 = 0.030069;
M  = M1*(1-x) + M2*x;

% other calculated properties
s_T    = cp/T + mu_T^2/mu_x;
Xi     = (ALPHA*rho+KAPPA^2/mu_x)*rho*pc/(rhoc)^2;
kB     = 1.38064852e-23;          % [J/K]
T_x    = - mu_x/mu_T;
zeta_0 = 0.150*(1-x) + 0.190*x;   % [nm]
% Olchowy, Sengers. International Journal of Thermophysics, 1989.
Gamma0 = 0.0520*(1-x)+ 0.0563*x;  
qD_1   = 0.5056*(1-x)+ 0.5827*x;  % [nm]
[alpha_b,beta_b,gamma_b,eta] = get_bg(x,T,Tc,rho,rhoc,p,pc,M);
       
% Eq.(17) in Kiselev, Huber. Fluid phase equilibria, 1998, 142.
zeta_OZ  = zeta_0*sqrt(Xi/Gamma0);
% Eq.(16)
zeta_hat = zeta_OZ*exp(-qD_1/zeta_OZ);
% Eq.(15)
z = zeta_hat/qD_1/10;
phi0 = 3/4/z^2*(1+z^2+(z^3-1/z)*atan(z))/(1+z^2);
% Eq.(14)
y1  = kB*T^2*rho/6/pi/eta/zeta_hat*1e9/gamma_b*s_T;
% Eq.(13)
y0  = kB*T*rho/6/pi/eta/zeta_hat*1e9/alpha_b/mu_x;
% Eq.(11)
yD  = 6*pi*eta^2/kB/T/rho*qD_1*1e-9/(phi0+1/y0);
% Eq.(12)
y1D = 6*pi*eta^2/kB/T/rho*qD_1*1e-9/(phi0+1/y1);
% Eq.(9)
OMEGA_alpha = 2/pi*(atan(zeta_hat/qD_1) - 1/sqrt(1+yD*zeta_hat/qD_1)*atan((zeta_hat/qD_1)/sqrt(1+yD*zeta_hat/qD_1)));
% Eq.(10)
OMEGA       = 2/pi*(atan(zeta_hat/qD_1) - 1/sqrt(1+y1D*zeta_hat/qD_1)*atan((zeta_hat/qD_1)/sqrt(1+y1D*zeta_hat/qD_1)));
% Eq.(6)
alpha = kB*T*rho/6/pi/eta/zeta_hat*1e9/mu_x*OMEGA_alpha + alpha_b;
% Eq.(7)
beta  = kB*T*rho/6/pi/eta/zeta_hat*1e9/T_x*OMEGA_alpha + beta_b;
% Eq.(23)
D = alpha/rho*mu_x;
% Eq.(26)
kT = T/rho*(alpha*mu_T+beta)/D;
% Eq.(20)
y = kB*T*rho/6/pi/eta/zeta_hat*1e9/alpha_b/mu_x*OMEGA_alpha;
y_star = beta_b/mu_T/alpha_b;
% Eq.(21)
Q = (y*(1+2*y_star)-y_star^2)/(1+y);
% Eq.(19)
lamda = kB*T*rho/6/pi/eta/zeta_hat*1e9*cp*OMEGA + alpha_b*mu_T^2*T*Q + gamma_b;

end