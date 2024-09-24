function [ylist]=mSciTrackCurve(f,y0,ytan,stepsize,nmax,tol,maxit,jach)
%f - function for bif diagram, one varying parameter must be an argument
%df - jacobian of said function at y0
%y0 - initial guess for first point on bifurcation curve 
%ytan - initial search direction for more equilibria (unit length)
%stepsize - length of each search vector
%nmax - number of equilibria to put on bif curve
%the rest: see mSciSolve

%Initialise ylist

ylength = length(y0); %our "n+1" from the worksheet
ylist = zeros(ylength,nmax);
%Find our first equilibrium:

fparamfix = @(x)f([x;y0(ylength)]);
x = mSciSolve(fparamfix,y0(1:(ylength-1)),tol,maxit,jach); 
disp('check')
y = [x;y0(ylength)]; %Append parameter value to the found equilibrium
ylist(:,1) = y;

for i = 2:nmax
    ypred = y + stepsize*ytan; %Initial guess for the next equilibrium
    ftosolve = @(y)[f(y);dot(ytan,y-ypred)]; %The system we aim to solve
    y = mSciSolve(ftosolve,y,tol,maxit,jach); %Solve for our next eq
    J = mSciJacobian(f,y,jach);
    ylist(:,i) = y;
    %Now to find the new search direction, ytan
    %We need the jacobian of f with fixed param, this is J from above
    %without the bottom row or rightmost column
    %size(J)
    %J = J(1:ylength-1,1:ylength-1);
    %now solve system in worksheet:
    %size([J;ytan'])
    %size([zeros(ylength-1,1);1]) 
    z = [J;ytan']\[zeros(ylength-1,1);1];
    ytan = (z/(max(abs(z))))*sign(z'*ytan);
end
    



