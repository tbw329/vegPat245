%x(2:4): klausmeierPO variables; x(1): time to run system for, x(5): P 
N = 1000;

%Other Params (from Gandhi et al)
D = 10; %Seed dispersal rate
L = 4; %Evap Rate
M = 1.8; %Mortality Rate
J = 0.003; %Water Use Efficiency
V = 63; %Advection Rate
R = 100; %Transpiration Rate

%Hardcoding this for finding eigenvalues at Hopf point for now
k = 0.244;

%Params for curve tracking alg
stepsize = 0.002;
nmax = 1000;
%These params are for rootfinding
tol = 1e-2;
maxit = 1000;
%This is for the Jacobian
jach = 1e-5;

%hard-coding in the approx location of the hopf point [b,w,beta,P]
PHopf = 272.6; %(From analysing busse balloon)
Beq = @(P) 0.5*(P*J/M + sqrt((P*J/M)^2-4*L/R));
bHopf = Beq(PHopf);
wHopf = M/(J*R*bHopf);

hopfPoint = [wHopf;bHopf;0;PHopf];

flowMap = @(x)(mSciIVP(@klausmeierPO,[wHopf;x(2);x(3)],x(4),[0 x(1)],N)); 
systemToSolve = @(x)(flowMap(x) - [wHopf;x(2);x(3)]);


%Jacobian eigenvalues at hopf point

A = zeros(2,2);

A(1,1) = V*1i*k - R*bHopf^2 - L;
A(1,2) = -2*R*wHopf*bHopf;
A(2,1) = J*R*bHopf^2;
A(2,2) = -D*k^2 + 2*J*R*wHopf*bHopf - M;

eval = eig(A);
eval = eval(2);

%Take imaginary part of eigenvalues to find frequency of local limit cycles
freq=abs(max(imag(eval)));
%Initial guess for period based on frequency of local limit cycles
T=2*pi/freq;

%Perturb the free variables at the Hopf point by a small amount alpha =
%10e-4

alpha = 10e-4;

hopfPointPerturbed = hopfPoint + [0;alpha;alpha;alpha];

y0 = [T;hopfPointPerturbed(2:4)];

%hardcoding initial ytan:
ytan = [0;1;0;0];


ylist = mSciTrackCurve(systemToSolve,y0,ytan,stepsize,nmax,tol,maxit,jach);

plot(ylist(5),ylist(2))
