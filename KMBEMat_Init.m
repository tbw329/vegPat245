%KlausMeier Backward Euler Matrix Initialisation
%Again, this is being written so the matrices can be calculated outside of f_oneDPDE.m
%dt, dx are the timestep and gridstep sizes, 
%V and D are parameters, 
% xmax is the total number of gridsteps (domain length is dx*xmax)

function [lambdaMatInv, muMatInv] = KMBEMat_Init(dt,dx,V,D,xmax)

lambda = dt*V/dx;
mu = D*dt/(dx^2);

lambdaMat = diag(repmat(1+lambda,xmax,1)) + diag(repmat(-lambda,xmax-1,1),1);
lambdaMat(xmax,1) = -lambda;
lambdaMatInv = lambdaMat^-1;

muMat = diag(repmat(1+2*mu,xmax,1)) + diag(repmat(-mu,xmax-1,1),1) + diag(repmat(-mu,xmax-1,1),-1);
muMat(xmax,1) = -mu;
muMat(1,xmax) = -mu;
muMatInv = muMat^-1;