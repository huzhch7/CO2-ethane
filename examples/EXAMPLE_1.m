clc
clear

%% input parameters
p_ref = 8.0e6;    % Pressure,      [Pa]
T_ref = 316.39;   % Temperature,   [K]
c_ref = 0.6721;   % Mass fraction of C2H6, [-]

%% obtain properties
prop = CO2_C2H6(p_ref,T_ref,c_ref);

%% show results
result{1,1} = '**Reference state**';
result{1,2} = '********';
result{1,3} = '********';
result{1,4} = '********';

result{2,1} = 'T';
result{2,2} = T_ref;
result{2,3} = 'K';
result{2,4} = 'Temperature';

result{3,1} = 'p';
result{3,2} = p_ref/1e6;
result{3,3} = 'MPa';
result{3,4} = 'Pressure';

result{4,1} = 'c';
result{4,2} = c_ref;
result{4,3} = '-';
result{4,4} = 'Mass fraction C2H6';

result{5,1} = '**Critical Parameters**';
result{5,2} = '********';
result{5,3} = '********';
result{5,4} = '********';

result{6,1} = 'Tc';
result{6,2} = prop.Tc;
result{6,3} = 'K';
result{6,4} = 'Critical temperature';

result{7,1} = 'pc';
result{7,2} = prop.pc/1000000;
result{7,3} = 'MPa';
result{7,4} = 'Critical pressure';

result{8,1} = 'rhoc';
result{8,2} = prop.rhoc;
result{8,3} = 'kg/m3';
result{8,4} = 'Critical density';

result{9,1} = '**Thermodynamic Properties**';
result{9,2} = '********';
result{9,3} = '********';
result{9,4} = '********';

result{10,1} = 'rho';
result{10,2} = prop.rho;
result{10,3} = 'kg/m3';
result{10,4} = 'Density';

result{11,1} = 'alpha';
result{11,2} = prop.alpha;
result{11,3} = 'Pa-1';
result{11,4} = 'Isothermal Compressibility';

result{12,1} = 'beta';
result{12,2} = prop.beta;
result{12,3} = 'K-1';
result{12,4} = 'Thermal Expansion coefficient';

result{13,1} = 'kappa';
result{13,2} = prop.kappa;
result{13,3} = '-';
result{13,4} = 'Concentration contraction coefficient';

result{14,1} = 'cs';
result{14,2} = prop.cs;
result{14,3} = 'kg/J';
result{14,4} = 'Concentration suspecbility';


result{15,1} = 'H';
result{15,2} = prop.H;
result{15,3} = 'J/kg';
result{15,4} = 'Partial enthalpy';

result{16,1} = 'M';
result{16,2} = prop.M;
result{16,3} = 'kg/mol';
result{16,4} = 'Molar mass';

result{17,1} = 'cp';
result{17,2} = prop.cp;
result{17,3} = 'J/(kg·K)';
result{17,4} = 'Specific heat';

result{18,1} = '**Transport Properties**';
result{18,2} = '********';
result{18,3} = '********';
result{18,4} = '********';

result{19,1} = 'eta';
result{19,2} = prop.eta;
result{19,3} = 'Pa·s';
result{19,4} = 'Viscosity';

result{20,1} = 'lambda';
result{20,2} = prop.lam;
result{20,3} = 'W/(m·K)';
result{20,4} = 'Thermal conductivity';

result{21,1} = 'D';
result{21,2} = prop.D;
result{21,3} = 'm^2/s';
result{21,4} = 'Diffusion coefficient';

result{22,1} = 'kT';
result{22,2} = prop.kT;
result{22,3} = '-';
result{22,4} = 'Thermal diffusion factor';


% output
result
