%Set up 

%model (As in Gandhi et al)
D = 1; %Seed dispersal rate 
L = 4; %Evap Rate
M = 1.8; %Mortality Rate
J = 0.003; %Water Use Efficiency
V = 63; %Advection Rate
R = 100; %Transpiration Rate
P = 250; %Precip rate 

%Finite difference parameters
tmax = 10000; %number of timesteps to run
xmax = 100; %number of gridpoints in space
dt = 0.01; %Timestep size
dx = 1; %Grid resolution
xax = linspace(0,(xmax-1)*dx,xmax);

%Preallocating size of W,B
W = zeros(tmax,xmax);
B=W;

% %DBCs
% Wx0 = 0; 
% Bx0 = 0.5; 

%IC  (The exact values of these are chosen from plot in Gandhi et al), B=W for ICs. 

    %Impulse ICs
%     W(1,1) = 2;
%     B(1,1) = 2;

for j = 1:xmax

%         %Uniform Noisy ICs
%         B(1,j) = 0.3 + (0.3*rand-0.15); % Initial vegetation density
%         W(1,j) = 0.3 + (0.3*rand-0.15); % Initial water concentration

% %             Periodic ICs
%             B(1,j) = 0.3 + 0.15*(sin(2*pi*j/xmax)); % Initial vegetation density
%             W(1,j) = 0.3 + 0.15*(sin(2*pi*j/xmax)); % Initial water concentration

        %Step ICs
        if j<xmax/2
            B(1,j) = 0.3;
            W(1,j) = 0.3;
        else
            B(1,j) = 0.45;
            W(1,j) = 0.45;
        end

end
    
    figure(1)
    clf
    fig1 = gca;
    pl = plot(fig1,xax,B(1,:));
    
    ylim([0 2])
    drawnow

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


for i = 2:tmax
    Wrow = W(i-1,:);
    Brow = B(i-1,:);
    W(i,:) = (lambdaMatInv*(Wrow + dt*(L*Wrow - R*Wrow.*Brow.^2 + P))')';
    B(i,:) = (muMatInv*(Brow + dt*(J*R*Wrow.*Brow.^2 - M*Brow))')';

if mod(i,10) == 0
    pl.YData = B(i,:); 
    drawnow
    i
end
end

        
