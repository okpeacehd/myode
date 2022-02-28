%% Initial conditions
clear
clc
tspan = [0,45000];       % time span 
C0 = [1,0,0,0,0,0,0,0,0,1,0];  % initial concentrations for main
%% ODE Solver
[t,C] = ode45(@InternalODE,tspan,C0);
%% Rh

df=1.8;
mwM=150e3;
A = 5.4./(mwM.^(1/df));
n = size(C,2)-2;

n1 = (C(:,1).*(7.5).*((mwM).^2))+(C(:,2).*(5.4).*((mwM).^2));
d1 = (C(:,1).*mwM.^2)+(C(:,2).*mwM.^2);

denominator = d1 + sum(((C(:,(3:n))).*((((3:n)-1).*mwM).^2)),2);
numerator = n1 + sum((C(:,(3:n)).*A.*(((3:n)-1).*mwM).^(1/df).*((((3:n)-1).*mwM).^2)),2);
Rh = numerator./denominator;
     

figure(2)
plot(t,Rh);
%% Function
function dCdt = InternalODE(~,C)
%% Rates
% unfolding rate
kuf = 0.000751334761923340;
% Equilibrium constant
Keq = 0.793500293693322;
% folding rate
kf  = kuf./Keq;
%% Diffusion Limited Cluster Aggregation kernel
df = 1.8;                             % 1.8 - 2.4 fractal dimension
L  = 1-(1/df);                        % fractal aggregate scales
kb = 1.38*10^-23;                     % Boltzman constant (1.38*10^-23)
T  = 300;                             % Absolute temperature (K)
n  = 0.00089;                         % Viscosity (Pa s or kg m^-1 s^-2)
W  = 7.98554851370698e-15;            % Fitting parameter
                                      
ks = (8/3) .* (kb.*T)./n;             % Aggregation rate constant
B  = @(i, j) 1/4 .* ((i^(1/df))  ...  % as a function of mass i,j 
    + j^(1/df)).* (1/(i^(1/df))  ...
    +1/(j^(1/df)));
P = @(i, j) (i*j).^L;                 % Product kernel (additional Factor)
k = @(i, j) (ks./W) .* B(i,j) ...     % Aggregation rate via DLCA
    .* P(i,j);
%% Rate of Aggregation 
% Monomer unfolded
unfolded = -kf.*C(1)+kuf.*C(2);
% Monomer folded
folded =   kf.*C(1)-kuf.*C(2)-2.*(C(2)^2).*k(1,1)-(C(2).*C(3)).*k(1,2)-(C(2).*C(4))...
    .*k(1,3)-(C(2).*C(5)).*k(1,4)-(C(2).*C(6)).*k(1,5)-(C(2).*C(7)).*k(1,6)-(C(2).*C(8)).*k(1,7);
% Dimer    (2)
dimer =    (C(2)^2).*k(1,1)-(C(2).*C(3)).*k(1,2)-2.*(C(3).*C(3)).*k(2,2)-(C(4).*C(3))...
    .*k(3,2)-(C(5).*C(3)).*k(4,2)-(C(6).*C(3)).*k(5,2)-(C(7).*C(3)).*k(6,2);
% Trimer   (3)
trimer =   (C(2).*C(3)).*k(1,2)-(C(2).*C(4)).*k(1,3)-(C(3).*C(4)).*k(2,3)-2.*...
    (C(4).*C(4)).*k(3,3)-(C(5).*C(4)).*k(4,3)-(C(6).*C(4)).*k(5,3);
% Tetramer (4)
tetramer = (C(2).*C(4)).*k(1,3)+(C(3)^2).*k(2,2)-(C(2)*C(5).*k(1,4))-(C(3).*C(5).*k(2,4))...
    -(C(4).*C(5).*k(3,4))-2.*(C(5).*C(5).*k(4,4));
% pentamer (5)
pentamer = (C(2).*C(5)).*k(1,4)+(C(3).*C(4).*k(2,3))-(C(2).*C(6).*k(1,5))-(C(3).*...
    C(6).*k(2,5))-(C(4).*C(6).*k(3,5));
% hexamer  (6)
hexamer =  (C(2).*C(6)).*k(1,5)+(C(3).*C(5).*k(2,4))+(C(4).*C(4).*k(3,3))-(C(2).*...
    C(7).*k(1,6))-(C(3).*C(7).*k(2,6));
% heptamer (7)
heptamer = (C(2).*C(7)).*k(1,6)+(C(3).*C(6).*k(2,5))+(C(4).*C(5).*k(3,4))-(C(2).*C(8).*k(1,7));
% octamer  (8)
octamer =  (C(2).*C(8)).*k(1,7) + (C(3).*C(7).*k(2,6))+(C(4).*C(6).*k(3,5))+(C(5).*C(5).*k(4,4));

total = unfolded + folded + 2*dimer + 3*trimer + 4*tetramer + 5*pentamer +...
    6*hexamer + 7*heptamer + 8*octamer;

aggregated = 2*dimer + 3*trimer + 4*tetramer + 5*pentamer +...
    6*hexamer + 7*heptamer + 8*octamer;

dCdt=[unfolded; folded; dimer; trimer; tetramer; pentamer; hexamer; ...
    heptamer; octamer; total; aggregated];
end
