%% CALCULATE PROPERTIES AND PLOT
% code for Fig. 3 in Zhanchao Hu paper

clc
clear

%% set pressure and temperature
% the range of p
p_range = ones(4,200);
p_range(1,:) = p_range(1,:)*5.6e6;
p_range(2,:) = p_range(2,:)*6e6;
p_range(3,:) = p_range(3,:)*7e6;
p_range(4,:) = p_range(4,:)*8e6;
% the range of T
T_range = linspace(287,338,200);
T_range = repmat(T_range,4,1);

%% calculate properties
for i = 1:4
    for j = 1:200 
        % output process
        disp(['i = ',num2str(i),' to 4 ; j = ',num2str(j),' to 200 ;']);
        % c=0.6721
        prop = CO2_C2H6(p_range(i,j),T_range(i,j),0.6721); 
        % store results
        rho(i,j)  = prop.rho;
        alpha(i,j)  = prop.alpha;
        beta(i,j)  = prop.beta;
        kappa(i,j)  = prop.kappa;
        cp(i,j)  = prop.cp;
        H(i,j)  = prop.H;
        cs(i,j)  = prop.cs; 
        eta(i,j)  = prop.eta;
        kT(i,j)  = prop.kT;
        lam(i,j) = prop.lam;
        D(i,j) = prop.D;
    end
end

%% rho
figure('units','inches','position',[5 5 3.3 1.5])
plot(T_range(1,:),rho(1,:),'b-','linewidth',1.2)
hold on
plot(T_range(2,:),rho(2,:),'r--','linewidth',1.2)
plot(T_range(3,:),rho(3,:),'-.','linewidth',1.2,'color',[46 139 87]/255);
plot(T_range(4,:),rho(4,:),'k:','linewidth',1.2)
axis([280 340 0 500])
set(gca,'ticklabelinterpreter','latex')
xlabel('$T/({\rm K})$','interpreter','latex')
ylabel('$\rho/({\rm kg/m^3})$','interpreter','latex')
h = legend('$p=5.6$~MPa','$p=6$~MPa','$p=7$~MPa','$p=8$~MPa'); 
set(h,'interpreter','latex')

%% alpha
figure('units','inches','position',[5 5 3.3 1.5])
plot(T_range(1,:),alpha(1,:),'b-','linewidth',1.2)
hold on
plot(T_range(2,:),alpha(2,:),'r--','linewidth',1.2)
plot(T_range(3,:),alpha(3,:),'-.','linewidth',1.2,'color',[46 139 87]/255);
plot(T_range(4,:),alpha(4,:),'k:','linewidth',1.2)
axis([280 340 0 0.25e-5])
set(gca,'ticklabelinterpreter','latex')
xlabel('$T/({\rm K})$','interpreter','latex')
ylabel('$\alpha/({\rm Pa^{-1}})$','interpreter','latex')
h = legend('$p=5.6$~MPa','$p=6$~MPa','$p=7$~MPa','$p=8$~MPa'); 
set(h,'interpreter','latex')

%% beta
figure('units','inches','position',[5 5 3.3 1.5])
plot(T_range(1,:),beta(1,:),'b-','linewidth',1.2)
hold on
plot(T_range(2,:),beta(2,:),'r--','linewidth',1.2)
plot(T_range(3,:),beta(3,:),'-.','linewidth',1.2,'color',[46 139 87]/255);
plot(T_range(4,:),beta(4,:),'k:','linewidth',1.2)
axis([280 340 0 0.3])
set(gca,'ticklabelinterpreter','latex')
xlabel('$T/({\rm K})$','interpreter','latex')
ylabel('$\beta/({\rm K^{-1}})$','interpreter','latex')
h = legend('$p=5.6$~MPa','$p=6$~MPa','$p=7$~MPa','$p=8$~MPa'); 
set(h,'interpreter','latex')
 
%% kappa
figure('units','inches','position',[5 5 3.3 1.5])
plot(T_range(1,:),kappa(1,:),'b-','linewidth',1.2)
hold on
plot(T_range(2,:),kappa(2,:),'r--','linewidth',1.2)
plot(T_range(3,:),kappa(3,:),'-.','linewidth',1.2,'color',[46 139 87]/255);
plot(T_range(4,:),kappa(4,:),'k:','linewidth',1.2)
axis([280 340 0 10])
set(gca,'ticklabelinterpreter','latex')
xlabel('$T/({\rm K})$','interpreter','latex')
ylabel('$\kappa$','interpreter','latex')
h = legend('$p=5.6$~MPa','$p=6$~MPa','$p=7$~MPa','$p=8$~MPa'); 
set(h,'interpreter','latex')
 
%% kT
figure('units','inches','position',[5 5 3.3 1.5])
plot(T_range(1,:),kT(1,:),'b-','linewidth',1.2)
hold on
plot(T_range(2,:),kT(2,:),'r--','linewidth',1.2)
plot(T_range(3,:),kT(3,:),'-.','linewidth',1.2,'color',[46 139 87]/255);
plot(T_range(4,:),kT(4,:),'k:','linewidth',1.2)
axis([280 340 -1 13])
set(gca,'ticklabelinterpreter','latex')
xlabel('$T/({\rm K})$','interpreter','latex')
ylabel('$k_T$','interpreter','latex')
h = legend('$p=5.6$~MPa','$p=6$~MPa','$p=7$~MPa','$p=8$~MPa'); 
set(h,'interpreter','latex')

%% cp
figure('units','inches','position',[5 5 3.3 1.5])
plot(T_range(1,:),cp(1,:),'b-','linewidth',1.2)
hold on
plot(T_range(2,:),cp(2,:),'r--','linewidth',1.2)
plot(T_range(3,:),cp(3,:),'-.','linewidth',1.2,'color',[46 139 87]/255);
plot(T_range(4,:),cp(4,:),'k:','linewidth',1.2)
axis([280 340 1 4e4])
set(gca,'ticklabelinterpreter','latex')
xlabel('$T/({\rm K})$','interpreter','latex')
ylabel('$c_p/[{\rm J/(kg\cdot K)}]$','interpreter','latex')
h = legend('$p=5.6$~MPa','$p=6$~MPa','$p=7$~MPa','$p=8$~MPa'); 
set(h,'interpreter','latex')

%% cs
figure('units','inches','position',[5 5 3.3 1.5])
plot(T_range(1,:),cs(1,:),'b-','linewidth',1.2)
hold on
plot(T_range(2,:),cs(2,:),'r--','linewidth',1.2)
plot(T_range(3,:),cs(3,:),'-.','linewidth',1.2,'color',[46 139 87]/255);
plot(T_range(4,:),cs(4,:),'k:','linewidth',1.2)
axis([280 340 2e-6 1.2e-5])
set(gca,'ticklabelinterpreter','latex')
xlabel('$T/({\rm K})$','interpreter','latex')
ylabel('$c_s/({\rm kg/J})$','interpreter','latex')
h = legend('$p=5.6$~MPa','$p=6$~MPa','$p=7$~MPa','$p=8$~MPa'); 
set(h,'interpreter','latex')

%% eta
figure('units','inches','position',[5 5 3.3 1.5])
plot(T_range(1,:),eta(1,:),'b-','linewidth',1.2)
hold on
plot(T_range(2,:),eta(2,:),'r--','linewidth',1.2)
plot(T_range(3,:),eta(3,:),'-.','linewidth',1.2,'color',[46 139 87]/255);
plot(T_range(4,:),eta(4,:),'k:','linewidth',1.2)
axis([280 340 0 7e-5])
set(gca,'ticklabelinterpreter','latex')
xlabel('$T/({\rm K})$','interpreter','latex')
ylabel('$\eta/({\rm Pa\cdot s})$','interpreter','latex')
h = legend('$p=5.6$~MPa','$p=6$~MPa','$p=7$~MPa','$p=8$~MPa'); 
set(h,'interpreter','latex')
 

%% D
figure('units','inches','position',[5 5 3.3 1.5])
plot(T_range(1,:),D(1,:),'b-','linewidth',1.2)
hold on
plot(T_range(2,:),D(2,:),'r--','linewidth',1.2)
plot(T_range(3,:),D(3,:),'-.','linewidth',1.2,'color',[46 139 87]/255);
plot(T_range(4,:),D(4,:),'k:','linewidth',1.2)
axis([280 340 2e-8 1.4e-7])
set(gca,'ticklabelinterpreter','latex')
xlabel('$T/({\rm K})$','interpreter','latex')
ylabel('$D/({\rm m^2/s})$','interpreter','latex')
h = legend('$p=5.6$~MPa','$p=6$~MPa','$p=7$~MPa','$p=8$~MPa'); 
set(h,'interpreter','latex')

%% lambda
figure('units','inches','position',[5 5 3.3 1.5])
plot(T_range(1,:),lam(1,:),'b-','linewidth',1.2)
hold on
plot(T_range(2,:),lam(2,:),'r--','linewidth',1.2)
plot(T_range(3,:),lam(3,:),'-.','linewidth',1.2,'color',[46 139 87]/255);
plot(T_range(4,:),lam(4,:),'k:','linewidth',1.2)
axis([280 340 0 0.2])
set(gca,'ticklabelinterpreter','latex')
xlabel('$T/({\rm K})$','interpreter','latex')
ylabel('$\lambda/[{\rm W/(m\cdot K)}]$','interpreter','latex')
h = legend('$p=5.6$~MPa','$p=6$~MPa','$p=7$~MPa','$p=8$~MPa'); 
set(h,'interpreter','latex')

%% H
figure('units','inches','position',[5 5 3.3 1.5])
plot(T_range(1,:),H(1,:),'b-','linewidth',1.2)
hold on
plot(T_range(2,:),H(2,:),'r--','linewidth',1.2)
plot(T_range(3,:),H(3,:),'-.','linewidth',1.2,'color',[46 139 87]/255);
plot(T_range(4,:),H(4,:),'k:','linewidth',1.2)
axis([280 340 -1.3e6 0.2e6])
set(gca,'ticklabelinterpreter','latex')
xlabel('$T/({\rm K})$','interpreter','latex')
ylabel('$\bar H/({\rm J/kg})$','interpreter','latex')
h = legend('$p=5.6$~MPa','$p=6$~MPa','$p=7$~MPa','$p=8$~MPa'); 
set(h,'interpreter','latex')