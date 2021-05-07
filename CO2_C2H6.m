function [prop_new] = CO2_C2H6(p_value,T_value,c_value)

[dx,dT,dc] = deal(1e-5,1e-4,1e-5);

% mole fraction
M1 = 0.04401;  % molar mass of CO2, [kg/mol]
M2 = 0.030069; % molar mass of C2H6,[kg/mol]
% molar fraction of C2H6, [-]
x_value = c_value.*M1./(c_value.*M1+(1-c_value).*M2);

%% Mixture model of methane and ethane
HEOS1 = py.CoolProp.CoolProp.AbstractState('HEOS','Ethane&CarbonDioxide');
HEOS1.set_mole_fractions([x_value,1-x_value]);
HEOS1.specify_phase(py.CoolProp.CoolProp.iphase_supercritical);
HEOS1.update(py.CoolProp.CoolProp.PT_INPUTS, p_value, T_value);
% dx changes
HEOS2 = py.CoolProp.CoolProp.AbstractState('HEOS','Ethane&CarbonDioxide');
HEOS2.set_mole_fractions([x_value + dx,1 - x_value - dx]);
HEOS2.specify_phase(py.CoolProp.CoolProp.iphase_supercritical);
HEOS2.update(py.CoolProp.CoolProp.PT_INPUTS, p_value, T_value);
% dT changes
HEOS3 = py.CoolProp.CoolProp.AbstractState('HEOS','Ethane&CarbonDioxide');
HEOS3.set_mole_fractions([x_value,1-x_value]);
HEOS3.specify_phase(py.CoolProp.CoolProp.iphase_supercritical);
HEOS3.update(py.CoolProp.CoolProp.PT_INPUTS, p_value, T_value + dT);
% dc changes
prop.M = HEOS1.molar_mass();  % [kg/mol]
dx2 = prop.M^2/M1/M2*dc;
HEOS4 = py.CoolProp.CoolProp.AbstractState('HEOS','Ethane&CarbonDioxide');
HEOS4.set_mole_fractions([x_value + dx2,1 - x_value - dx2]);
HEOS4.specify_phase(py.CoolProp.CoolProp.iphase_supercritical);
HEOS4.update(py.CoolProp.CoolProp.PT_INPUTS, p_value, T_value);

%% Thermodynamics properties
prop.x = x_value;
prop.T = T_value;               % [K]
prop.p = p_value;               % [pa]
prop.rho = HEOS1.rhomolar();    % mol/m3
prop.rhomass = HEOS1.rhomass(); % kg/m3
prop.cv  = HEOS1.cvmolar();     % J/(mol K)
prop.cp  = HEOS1.cpmolar();     % J/(mol K)
prop.alpha = HEOS1.isothermal_compressibility();     %[Pa-1]       
prop.beta  = HEOS1.isobaric_expansion_coefficient(); %[K-1]
% kappa: dc base, used in flow governing equations
prop.kappa1 = (HEOS4.rhomolar()-HEOS1.rhomolar())/dc/HEOS1.rhomolar(); %[-]
% kappa: dx base, used in transport properties calculations
prop.kappa2 = (HEOS2.rhomolar()-HEOS1.rhomolar())/dx/HEOS1.rhomolar(); %[-]
prop.smolar_residual = HEOS1.smolar_residual();                  %[J/(mol·K)]

%% chemical potential, mass fraction base, used in flow governing equations
mu1 = HEOS1.chemical_potential(0)/M2-HEOS1.chemical_potential(1)/M1; % [J/kg]
mu3 = HEOS3.chemical_potential(0)/M2-HEOS3.chemical_potential(1)/M1; % [J/kg]
mu4 = HEOS4.chemical_potential(0)/M2-HEOS4.chemical_potential(1)/M1; % [J/kg]
prop.mu_c = (mu4-mu1)/dc;                % [J/kg]
prop.H    = mu1 - prop.T*((mu3-mu1)/dT); % [J/kg]

%% chemical potential: molar fraction base, used in transport properties calculations
mu1 = HEOS1.chemical_potential(0)-HEOS1.chemical_potential(1); % [J/mol]
mu2 = HEOS2.chemical_potential(0)-HEOS2.chemical_potential(1); % [J/mol]
mu3 = HEOS3.chemical_potential(0)-HEOS3.chemical_potential(1); % [J/mol]
prop.mu_x = (mu2-mu1)/dx;           % [J/mol]
prop.mu_T = (mu3-mu1)/dT;           % [J/mol K]

%% critical properties
CPs       = HEOS1.all_critical_points();
prop.Tc   = CPs(py.int(0)).T;           % [K]
prop.pc   = CPs(py.int(0)).p;           % [Pa]
prop.rhoc = CPs(py.int(0)).rhomolar;    % [mol/m3]

%% transport properties
[prop.eta,prop.lambda,prop.D,prop.kT] = tran_prop(prop.x,...
                                                  prop.T,...
                                                  prop.Tc,...
                                                  prop.rho*prop.M,...
                                                  prop.rhoc*prop.M,...
                                                  prop.p,...
                                                  prop.pc,...
                                                  prop.cp/prop.M,...    % [J/(kg·K)]
                                                  prop.mu_x/prop.M,...  % [J/kg]
                                                  prop.mu_T/prop.M,...  % [J/(kg·K)]
                                                  prop.alpha,...
                                                  prop.kappa2);

%% Thermodynamics properties, convert unit and store
% Molar mass, [kg/mol]
prop_new.M    = prop.M;
% Critical temperature, [K]
prop_new.Tc   = prop.Tc;
% Critical density, [kg/m^3]
prop_new.rhoc = prop.rhoc*prop.M;
% Critical pressure, [Pa]
prop_new.pc   = prop.pc;
% Density, [kg/m^3]
prop_new.rho  = prop.rho*prop.M;

% Isothermal compressibility, [Pa^{-1}]
prop_new.alpha = prop.alpha;
% Thermal expansion coefficient, [K^{-1}]
prop_new.beta  = prop.beta;
% concentration contraction coefficient, [-]
prop_new.kappa = prop.kappa1;
% isobaric specific heat, [J/mol] to [J/kg]
prop_new.cp = prop.cp/prop.M;
% residual entropy, [J/(mol·K)] to [J/(kg·K)]
prop_new.sr = prop.smolar_residual/prop.M;
% partial enthalpy, [J/kg]
prop_new.H  = prop.H;
% Concentration susceptibility, [-]
prop_new.cs = 1.0/prop.mu_c;

%% Transport properties, store
% Viscosity, [Pa·s]
prop_new.eta = prop.eta;
% Thermal diffusion factor, [-]
prop_new.kT  = prop.kT;
% Thermal conductivity, [W/m·K]
prop_new.lam = prop.lambda;
% Diffusion coefficient, [m^2/s]
prop_new.D   = prop.D;
                                             
end