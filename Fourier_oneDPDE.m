%Set up 
clear all

%model (As in Gandhi et al)
D = 10; %Seed dispersal rate - Something funky going on here
L = 4; %Evap Rate
M = 1.8; %Mortality Rate
J = 0.003; %Water Use Efficiency
V = 63; %Advection Rate
R = 100; %Transpiration Rate
P = 241; %Precip rate (in 240<P<272 range for pattern formation

%Finite difference parameters
tmax = 200000; %number of timesteps to run
xmax = 1000; %number of gridpoints in space
dt = 0.01; %Timestep size
dx = 0.1; %Grid resolution
xax = linspace(0,(xmax-1)*dx,xmax);
domlength = xmax*dx;

%Preallocating size of W,B
W = zeros(tmax,xmax);
B=W;
Bsum = B(:,1);

% %DBCs
% Wx0 = 0; 
% Bx0 = 0.5; 

%ICs (perturbations about equilibrium
Beq = 0.5*(P*J/M + sqrt((P*J/M)^2-4*L/R)); 
Weq = M/(J*R*Beq);

    %Impulse ICs
%     W(1,1) = 2;
%     B(1,1) = 2;

for j = 1:xmax

%         %Uniform Noisy ICs
%         B(1,j) = 0.3 + (0.3*rand-0.15); % Initial vegetation density
%         W(1,j) = 0.3 + (0.3*rand-0.15); % Initial water concentration

%             Periodic ICs
            
             B(1,j) = Beq + 0.1*(sin(2*pi*4*j/xmax)); % Initial vegetation density
             W(1,j) = Weq + 0.1*(sin(2*pi*4*j/xmax)); % Initial water concentration

%             B(1,j) = Beq + 0.01*(rand(1)-1/2); % Initial vegetation density
%             W(1,j) = Weq + 0.01*(rand(1)-1/2); % Initial water concentration
            

%         %Step ICs
%         if j<xmax/2
%             B(1,j) = 0.3;
%             W(1,j) = 0.3;
%         else
%             B(1,j) = 0.45;
%             W(1,j) = 0.45;
%         end

end

%Keeping note of total biomass
Bsum(1) = sum(B(1,:))/xmax;

    figure(1)
    clf
    fig1 = gca;
    xlabel('Position [m]')
    ylabel('Biomass Density [kg/m^{2}]')
    pl1 = plot(fig1,xax,B(1,:),'LineWidth',3);
    figure(3)
    clf
    fig3 = gca;
    xlabel('Fourier Mode')
    ylabel('Fourier Coeff')
    pl3 = plot(fig3,B(1,:),'.','MarkerSize',10);

    
    %ylim([0 2])
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
    %P = 270-i/1000; %linearly decreasing P option
    Wrow = W(i-1,:);
    Brow = B(i-1,:);
    W(i,:) = (lambdaMatInv*(Wrow + dt*(-L*Wrow - R*Wrow.*Brow.^2 + P))')';
    B(i,:) = (muMatInv*(Brow + dt*(J*R*Wrow.*Brow.^2 - M*Brow))')';

if mod(i,100) == 0
    pl1.YData = B(i,:); 
    drawnow
    [~,sines,cosines,~] = fourierFinder(B(i,:),domlength,20);
    pl3.YData = sines+cosines;
    ylim([-1 1]);
    drawnow
    i
end
if mod(i,2000) == 0
    P = P-1;
end
Bsum(i) = sum(B(i,:))/xmax;
end

xlabel('Position[m]','FontSize',24)
ylabel('Biomass Density [kg/m^{2}]','FontSize',24)
ax=gca;
ax.FontSize = 24;
xlabel('Position [m]')
ylim([0 1])

%Total biomass over time plot, interesting
figure(2)
fig2 = gca;
pl2 = plot(fig2,1:tmax,Bsum);
ylabel('Total Biomass')
xlabel('Time')
title('Total Biomass Over Time')



        
