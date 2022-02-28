%% Initial conditions
clear
clc
tspan = [0,45000];       % time span 
C0 = [1,0,0,0,0,0,0,0,0,1,0];  % initial concentrations for main

%% ODE Solver
[t,C] = ode45(@InternalODE,tspan,C0);

for Keq = 0.1:0.1:1
[t2,D] = ode45(@(t,C)Lines(t,C,Keq),tspan,[1,0]);
figure(4);
hold on;
D=D(:,[1,2]);
plot(t2,D);
end
%% rH calculation

df=1.8; % Fractal Dimension
% Defined Solved Concentration Matrices 
unfolded = C(:,1);          
folded = C(:,2);
dim =  C(:,3);
tri =  C(:,4);
tet =  C(:,5);
pent = C(:,6);
hex =  C(:,7);
hep =  C(:,8);
oct =  C(:,9); %8

%Rh species equations
RhMon = 7.5;  %nm
RhMon2 = 5.4;
mwM= 150e3;
A = RhMon./(mwM.^(1/df));
RhDim =  A.*(2*mwM).^(1/df);
RhTri =  A.*(3*mwM).^(1/df);
RhTet =  A.*(4*mwM).^(1/df);
RhPent = A.*(5*mwM).^(1/df);
RhHex =  A.*(6*mwM).^(1/df);
RhHep =  A.*(7*mwM).^(1/df);
RhOct =  A.*(8*mwM).^(1/df);

%total = unfolded + folded + 2*dim + 3*tri + 4*tet + 5*pent +...
    %6*hex + 7*hep + 8*oct;

%aggregated = 2*dim + 3*tri + 4*tet + 5*pent +...
    %6*hex + 7*hep + 8*oct;

Rh = (((unfolded).*(RhMon).*((mwM).^2))+((folded).*(RhMon2).*((mwM).^2))+((dim).*(RhDim).*((mwM*2).^2))+((tri).*...
   (RhTri).*(mwM*3).^2)+((tet).*(RhTet).*(mwM*4).^2)+((pent).*(RhPent).*(mwM*5).^2)+((hex).*...
   (RhHex).*(mwM*6).^2)+((hep).*(RhHep).*(mwM*7).^2)+((oct).*(RhOct).*(mwM*8).^2))./...
   (((unfolded).*(mwM).^2)+((folded).*((mwM).^2))+((dim).*(mwM*2).^2)+((tri).*...
   (mwM*3).^2)+((tet).*(mwM*4).^2)+((pent).*(mwM*5).^2)+((hex).*...
   (mwM*6).^2)+((hep).*(mwM*7).^2)+((oct).*(mwM*8).^2));


%% plot 
figure(1);
t_vec = 0:1200:6000;
C_des = [0 0 1; 0.139059662 0.46212413 0.3988162; 0.546943487, 0.270364756, 0.182691757; 0.715038793, 0.193128843, 0.091832364; 0.841549484, 0.088989431, 0.069461085; 0.870984328, 0.116516682, 0.01249899];
B = C(:,[11, 1, 2]);
plot(t,B,"DisplayName","Result")
hold on
scatter(t_vec,C_des,'bo', "DisplayName", "Desired");

figure(2);
grid on; 
plot(t,C);
xlabel("Time");
ylabel("concentration");
legend(["Folded","Unfolded","Dimer","Trimer","Tetramer","Pentamer","Hexamer","Heptamer","Octamer","Total","Aggregated"]);

figure(3);
t_rh = [294;1339;2409;3476;4533;5591;6667;7753;8841;9938;11072;12182;13303;...
14431;15547;16682;17808;18950;20103;21261;22427;23587;24751;25922;27099;...
28276;29424;30579;31749;32923;34088;35264;36458;37659;38859;40063;41301;...
42533];
C_rH = [7.345145473;8.492829616;12.75222956;14.65870449;16.15314224;...
17.00860117;17.93371182;18.29812154;18.81330567;19.79549585;20.40703612;...
20.28250728;20.43511705;21.25386888;21.13328606;22.02857219;21.74014617;...
22.24101453;22.19476358;22.23872033;22.60615839;22.98827928;22.96212548;...
23.36452704;22.86007974;23.63726079;24.15116016;23.47437304;24.43169416;...
23.38664307;24.16850427;23.97377311;24.80656541;24.56631745;25.50941458;...
24.86713212;24.98909146;25.11820867];

plot(t,Rh);
hold on
scatter(t_rh,C_rH,'bo', "DisplayName", "Desired");
xlabel("Time");
ylabel("rH");

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

function dCdt = Lines(~,C,Keq)
%% Rates
% unfolding rate
kuf = 1;
% Equilibrium constant
% folding rate
kf  = kuf./Keq;
%% Unfolding and folding
% Monomer unfolded
unfolded = -kf.*C(1)+kuf.*C(2);
% Monomer folded
folded =   kf.*C(1)-kuf.*C(2);

dCdt=[unfolded; folded];
end
