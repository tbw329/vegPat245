function [df] = klausmeierPO(x,P) %x = [w;b;beta;P]

%Other Params (from Gandhi et al)
D = 10; %Seed dispersal rate 
L = 4; %Evap Rate
M = 1.8; %Mortality Rate
J = 0.003; %Water Use Efficiency
V = 63; %Advection Rate
R = 100; %Transpiration Rate

%At the moment, fix c too
c = 10;

db = x(3);
dw = (1/(V+c))*(L*x(1)+R*x(1)*(x(2)^2) - P);
dbeta = (1/D)*(M*x(2)-J*R*x(1)*(x(2)^2)-c*x(3));
df = [dw;db;dbeta];