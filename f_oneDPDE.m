%Set up 
function [Bavg,Btmax,Wtmax] = f_oneDPDE(P,B0,W0)

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

%Preallocating size of W,B
W = zeros(tmax,xmax);
B=W;

%Initial Conditions
B(1,:) = B0;
W(1,:) = W0;
            
%Initialise the matrices

lambda = dt*V/dx;
mu = D*dt/(dx^2);

lambdaMat = diag(repmat(1+lambda,xmax,1)) + diag(repmat(-lambda,xmax-1,1),1);
lambdaMat(xmax,1) = -lambda;
lambdaMatInv = lambdaMat^-1;

muMat = diag(repmat(1+2*mu,xmax,1)) + diag(repmat(-mu,xmax-1,1),1) + diag(repmat(-mu,xmax-1,1),-1);
muMat(xmax,1) = -mu;
muMat(1,xmax) = -mu;
muMatInv = muMat^-1;

%Time Stepping

for i = 2:tmax
    Wrow = W(i-1,:);
    Brow = B(i-1,:);
    W(i,:) = (lambdaMatInv*(Wrow + dt*(-L*Wrow - R*Wrow.*Brow.^2 + P))')';
    B(i,:) = (muMatInv*(Brow + dt*(J*R*Wrow.*Brow.^2 - M*Brow))')';
end

Btmax = B(tmax,:);
Wtmax = W(tmax,:);
Bavg = sum(Btmax)/xmax;