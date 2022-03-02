%% Initial conditions
clear
clc
tspan = [0,45000];       % time span 
C0 = [1,0,0,0,0,0,0,0,0,1,0];  % initial concentrations for main
%% ODE Solver
[t,C] = ode45(@InternalODE,tspan,C0);
plot(t,C);
%% Rh
%constants
df=1.8;
mwM=150e3;
A = 5.4./(mwM.^(1/df));
n = size(C,2)-2;
% initials
n1 = (C(:,1).*(7.5).*((mwM).^2))+(C(:,2).*(5.4).*((mwM).^2));
d1 = (C(:,1).*mwM.^2)+(C(:,2).*mwM.^2);
%vectorised loop
denominator = d1 + sum(((C(:,(3:n))).*((((3:n)-1).*mwM).^2)),2);
numerator = n1 + sum((C(:,(3:n)).*A.*((((3:n)-1).*mwM).^(1/df)).*((((3:n)-1).*mwM).^2)),2);
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
x=9;
for i = 2:x
    for j = 2:x
        if i==j
            kArray(i,j) = 2.*k(i-1,j-1).*C(i).*C(j);
            bArray(i,j) = k(i-1,j-1).*C(i).*C(j);
        else 
            kArray(i,j) = k(i-1,j-1).*C(i).*C(j);
            bArray(i,j) = k(i-1,j-1).*C(i).*C(j);
        end
    end
end
     
% Monomer unfolded
unfolded = -kf.*C(1)+kuf.*C(2);
% Monomer folded
folded =  kf.*C(1)-kuf.*C(2)-sum(kArray(2,2:x-1));

% Dimer    (2)
dimer =    bArray(2,2)-sum(kArray(3,2:x-2));
% Trimer   (3)
trimer =   bArray(2,3)-sum(kArray(4,2:x-3));
% Tetramer (4)
tetramer = sum((bArray(2:3,4:-1:3)))-sum(kArray(5,2:x-4));
% pentamer (5)
pentamer = sum((bArray(2:3,5:-1:4)))-sum(kArray(6,2:x-5));
% hexamer  (6)
hexamer =  sum((bArray(2:4,6:-1:4)))-sum(kArray(7,2:x-6));
% heptamer (7)
heptamer = sum((bArray(2:4,7:-1:5)))-sum(kArray(8,2:x-7));
% octamer  (8)
octamer =  sum((bArray(2:5,8:-1:5)));

total = unfolded + folded + 2.*dimer + 3.*trimer + 4.*tetramer + 5.*pentamer +...
    6.*hexamer + 7.*heptamer + 8.*octamer;

aggregated = 2.*dimer + 3.*trimer + 4.*tetramer + 5.*pentamer +...
    6.*hexamer + 7.*heptamer + 8.*octamer;

dCdt=[unfolded; folded; dimer; trimer; tetramer; pentamer; hexamer; ...
    heptamer; octamer; total; aggregated];
end
