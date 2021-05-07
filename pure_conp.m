function [eta,lamda_0,lamda_ex] = pure_conp(T,rho,i)
% the thermodynamic properties of pure substances
% i = 1: CO2, i=2: C2H6

%% constants
epkB  = [195.2,215.7];
sigma = [0.3941,0.4443];
a = [-9.83026,-2.88240;274.922,-2197.33;12.0085,7.78198;-3900.84,23144.8;
      0.05188,0.54518;24.1829,-190.691;2391.80,16158.2];
a = [a(:,2),a(:,1)];
b = [16.4265,14.5845;10.7440,20.9392;-6.13841,-14.0533;0.33633,2.14884;
    1.78933,3.75428;-0.19234,-0.98048];
b = [b(:,2),b(:,1)];
M = [44.01,30.069];
rhoc = [10.625,6.8701];
Tc   = [304.1282,305.32];

%% viscosity, unit [Pa s]
t        = T/epkB(i);
OMEGA2_2 = 1.16145/t^0.14874 + 0.52487*exp(-0.77320*t) + 2.16178*exp(-2.43787*t);
eta_0    = 26.69167e-9*sqrt(M(i)*T)/OMEGA2_2/sigma(i)^2;
eta_ex   = 1e-7*exp( a(1,i)+a(2,i)/T )*( exp( (a(3,i) + a(4,i)/T^(3/2))*rho^0.1 + (rho/rhoc(i)-1)*rho^0.5*(a(5,i)+a(6,i)/T+a(7,i)/T^2) ) -1 );
eta      = eta_0 + eta_ex;

%% thermal conductivity, unit [W/m K]
switch i
    case 1
        lam_coef = [1.51874307e-2,2.80674040e-2,2.28564190e-2,-7.41624210e-3];
        SUM = 0;
        for ind = [1,2,3,4]
            SUM = SUM + lam_coef(ind)/(T/Tc(i))^(ind-1);
        end
        lamda_0 = sqrt(T/Tc(i))/SUM/1000;
    case 2
        cp0_R    = 4 + exp(-0.02*T/100) * (0.0036890096*(T/100)^(4-1) + ...
                                          -0.17196907*(T/100)^(4-2) + ...
                                           3.159226*(T/100)^(4-3)  + ...
                                          -8.0459942*(T/100)^(4-4) + ...
                                           7.4237673*(T/100)^(4-5) + ...
                                           0*(T/100)^(4-6) + ...
                                          -2.0724572*(T/100)^(4-7)); 
        G       = 0.444358*(T/264.70)^0 +  0.327867*(T/264.70)^1 + 0.1936835*(T/264.70)^2;   
        lamda_0 = 0.177586*sqrt(T/M(i))*cp0_R/0.43075^2/G/1000;
   end

lamda_ex = 1e-3*rho/rhoc(i)*( b(1,i) + b(2,i)*(rho/rhoc(i))^2 + (b(3,i)+b(4,i)*Tc(i)/T)*(rho/rhoc(i))^3 + (b(5,i)+b(6,i)*Tc(i)/T)*(rho/rhoc(i))^4 );

end