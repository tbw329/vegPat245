%I turned the backward Euler step into a function to decouple the time to
%run (tmax) from the numerical integration procsess

%muMatInv,lambdaMatInv are calculate in f_oneDPDE.m
%dt as well as the parameters (L,R,J,P,M) are specified in the arguments/code of f_oneDPDE.m 
function [W_out,B_out] = klausmeierBackwardEuler(W_in,B_in,L,R,J,P,M,dt,lambdaMatInv,muMatInv)

%Backward Euler Algorithm
W_out = (lambdaMatInv*(W_in + dt*(-L*W_in - R*W_in.*B_in.^2 + P))')';
B_out = (muMatInv*(B_in + dt*(J*R*W_in.*B_in.^2 - M*B_in))')';

end
