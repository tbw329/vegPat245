function df = mSciJacobian(f,x,h)
% f: Function to differentiate
% x: Point we're differentiating at
% h: Finite difference to use when numerically diffing

%Initialising df:
xlength = length(x);
flength = length(f(x));


df = zeros(flength,xlength);

for i = 1:xlength
    hvec = zeros(xlength,1); %Add finite difference
    hvec(i) = h;
    df(:,i) = (f(x+hvec)-f(x-hvec))/(2*h); %Leapfrog Formula
end

        