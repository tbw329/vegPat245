%This function aims to find the first nterms fourier coefficients of a wave B
%on a domain [0,domLength]
function [a0,sines,cosines Fseries] = fourierFinder(B,domLength,nterms)

L = length(B); %L is the number of points on the wave
dx = domLength/L; %step size for trapezium method
a0 = (1/(domLength))*(dx/2)*(B(1)+B(L)+2*sum(B(2:L-1))); %Use trapezium method to integrate

Fseries = a0; %start collating Fourier series

%Initialise sines and cosines
sines = zeros(1,nterms);
cosines = sines;
domainSteps = 0:dx:domLength-dx;

%Trapezium integration to find cos and sin coefficients
for i = 1:nterms 
    sinewave = sin(i*pi*domainSteps/(domLength/2));
    sinefunc = B.*sinewave; 
    cosinewave = cos(i*pi*domainSteps/(domLength/2));
    cosinefunc = B.*cosinewave;
    sines(i) = (1/(domLength/2))*(dx/2)*(sinefunc(1)+sinefunc(L)+2*sum(sinefunc(2:L-1)));
    cosines(i) = (1/(domLength/2))*(dx/2)*(cosinefunc(1)+cosinefunc(L)+2*sum(cosinefunc(2:L-1)));
    Fseries = Fseries + sines(i)*sinewave + cosines(i)*cosinewave; %Add new terms to fourier series
end


