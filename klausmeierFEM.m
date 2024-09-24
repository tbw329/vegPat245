%Set up PDE
klausmeierModelSize = 2;
klausmeier = createpde(klausmeierModelSize);

%Define Domain [0,10]x[0,10]
Lx = 10;
Ly = 10;
rect = [3 4 0 Lx Lx 0 0 0 Ly Ly]'; %10x10 Rectangle specification for decsg
g = decsg(rect, 'R1', ('R1')'); %Region called R1
geometryFromEdges(klausmeier, g); %Using decomposed shape w/ model to make geometry for PDE

%Params (As in Gandhi et al)
D = 10; %Seed dispersal rate 
L = 4; %Evap Rate
M = 1.8; %Mortality Rate
J = 0.003; %Water Use Efficiency
V = 63; %Advection Rate
R = 100; %Transpiration Rate
P = 230; %Precip rate (chosen arb)

%Coeffs for PDE
m = [1;1]; %(du/dt coeffs)
d = [0,0;0,0]; %(d2u/dt2 coeffs)
c = [0;0;0;D;0;D]; %(2nd spatial derivatives)
a = @(region,state) [L+R*(state.u(2,:)).^2;M-J*R*state.u(1,:).*state.u(2,:)]; %Terms with W or B coeffs
f = @(region,state) [P+V*state.ux(1,:);0];

specifyCoefficients(klausmeier, ...
    'm', m, ... % Mass coefficients for both equations
    'd', d, ... % No damping
    'c', c, ... % Diffusion coefficients
    'a', a, ... % Reaction terms
    'f', f); 

%Neumann BCs
applyBoundaryCondition(klausmeier, 'neumann', 'Edge', 1:klausmeier.Geometry.NumEdges, 'g', [0; 0], 'q', [0; 0]);

%ICs, random noise
W0 = @(location) [0.5 + 0.1 * rand(size(location.y)); % Initial vegetation density
  1.0 + 0.1 * rand(size(location.y))]; % Initial water concentration
setInitialConditions(klausmeier, W0,0);

%Mesh Generation
generateMesh(klausmeier, 'Hmax', 0.1); % Adjust 'Hmax' for mesh resolution

%Solving System
tlist = linspace(0, 1, 100); % Define time steps
result = solvepde(klausmeier, tlist);

%Extract Results
W = result.NodalSolution(:,:,1); % Solution for W
B = result.NodalSolution(:,:,2); % Solution forB

%Plot
figure;
pdeplot(klausmeier, 'XYData', B(:,end), 'Contour', 'on');
title('Vegetation Density at Final Time Step');


