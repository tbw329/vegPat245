%Set up 

%model (As in Gandhi et al)
D = 1; %Seed dispersal rate 
L = 4; %Evap Rate
M = 1.8; %Mortality Rate
J = 0.3; %Water Use Efficiency
V = 63; %Advection Rate
R = 100; %Transpiration Rate
P = 2.5; %Precip rate 

%Finite difference parameters
tmax = 10000; %number of timesteps to run
xmax = 100; %number of gridpoints in space
dt = 0.001; %Timestep size
dx = 1; %Grid resolution
xax = linspace(0,(xmax-1)*dx,xmax);

%Preallocating size of W,B
W = zeros(tmax,xmax);
B=W;

% %DBCs
% Wx0 = 0; 
% Bx0 = 0.5; 

%IC  (The exact values of these are chosen from plot in Gandhi et al), B=W for ICs. 

%     %Impulse ICs
%     W(1,1) = 2;
%     B(1,1) = 2;

for j = 1:xmax

        %Uniform Noisy ICs
        B(1,j) = 0.3 + (0.3*rand-0.15); % Initial vegetation density
        W(1,j) = 0.3 + (0.3*rand-0.15); % Initial water concentration

%         %Periodic ICs
%         B(1,j) = 0.3 + 0.15*(sin(2*pi*j/xmax)); % Initial vegetation density
%         W(1,j) = 0.3 + 0.15*(sin(2*pi*j/xmax)); % Initial water concentration

%         %Step ICs
%         if j<xmax/2
%             B(1,j) = 0.3;
%             W(1,j) = 0.3;
%         else
%             B(1,j) = 0.45;
%             W(1,j) = 0.45;
%         end

end

for i = 2:tmax
    i
     %Periodic BCs 
     W(i,1) = W(i-1,1)+dt*((V/(2*dx))*(W(i-1,2)-W(i-1,xmax))-R*W(i-1,1)*B(i-1,1)^2+P - L*W(i-1,1));
     B(i,1) = B(i-1,1)+dt*(D/(dx^2)*(B(i-1,2)-2*B(i-1,1)+B(i-1,xmax))-M*B(i-1,1)+J*R*W(i-1,1)*B(i-1,1)^2);
     W(i,xmax) = W(i-1,xmax)+dt*((V/(2*dx))*(W(i-1,1)-W(i-1,xmax-1))-R*W(i-1,xmax)*B(i-1,xmax)^2+P - L*W(i-1,xmax));
     B(i,xmax) = B(i-1,xmax)+dt*(D/(dx^2)*(B(i-1,1)-2*B(i-1,xmax)+B(i-1,xmax-1))-M*B(i-1,xmax)+J*R*W(i-1,xmax)*B(i-1,xmax)^2);
    for j = 2:xmax-1
        W(i,j) = W(i-1,j)+dt*((V/(2*dx))*(W(i-1,j+1)-W(i-1,j-1))-R*W(i-1,j)*B(i-1,j)^2+P - L*W(i-1,j));
        B(i,j) = B(i-1,j)+dt*(D/(dx^2)*(B(i-1,j+1)-2*B(i-1,j)+B(i-1,j-1))-M*B(i-1,j)+J*R*W(i-1,j)*B(i-1,j)^2);
    end
%     W(i,1) = W(i,2);
%     W(i,xmax) = W(i,xmax-1);
%     B(i,1) = B(i,2);
%     B(i,xmax) = B(i,xmax-1);
    plot(xax,B(i,:))
    ylim([0 2])
    drawnow
end

        
