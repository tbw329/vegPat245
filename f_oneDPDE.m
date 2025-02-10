%Set up 
function [Bavg,Btmax,Wtmax,xax] = f_oneDPDE(P,B0,W0)

%model (As in Gandhi et al, P now variable)
D = 10; %Seed dispersal rate - Something funky going on here
L = 4; %Evap Rate
M = 1.8; %Mortality Rate
J = 0.003; %Water Use Efficiency
V = 63; %Advection Rate
R = 100; %Transpiration Rate

xmax = length(B0); %number of gridpoints in space

%Other finite difference parameters, hard coding these at the moment
tmax = 2000; %number of timesteps to run
dt = 0.01; %Timestep size
dx = 0.1; %Grid resolution 
xax = linspace(0,(xmax-1)*dx,xmax);

%Preallocating size of W,B
W = zeros(tmax,xmax);
B=W;

%Initial Conditions
B(1,:) = B0;
W(1,:) = W0;
            
%Initialise the matrices


[lambdaMatInv, muMatInv] = KMBEMat_Init(dt,dx,V,D,xmax);

%Time Stepping

for i = 2:tmax
    [W(i,:),B(i,:)] = klausmeierBackwardEuler(W(i-1,:),B(i-1,:),L,R,J,P,M,dt,lambdaMatInv,muMatInv);
end

Btmax = B(tmax,:);
Wtmax = W(tmax,:);
Bavg = sum(Btmax)/xmax;