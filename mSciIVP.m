function [xend,t,xt]=mSciIVP(f,x0,P,tspan,N)

%Initialising
xt = zeros(length(f(x0,P)),N+1);
t = zeros(N+1,1);
t(1) = tspan(1);
xt(:,1) = x0;
h = (tspan(2)-tspan(1))/N;

for i = 1:N
    %RK4 Timestepping for autonomous system.
    k1 = f(xt(:,i),P);
    k2 = f(xt(:,i) + (1/2)*k1*h,P);
    k3 = f(xt(:,i) + (1/2)*k2*h,P);
    k4 = f(xt(:,i) + k3*h,P);

    xt(:,i+1) = xt(:,i) + h*((1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4);
    t(i + 1) = t(i) + h;
end

%Defining xend as last entry in numerical intergration vector.
xend = xt(:,N+1);