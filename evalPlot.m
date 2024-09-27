
clf
hold on
%Param values from Gandhi et al
D = 10; %Seed dispersal rate 
L = 4; %Evap Rate
M = 1.8; %Mortality Rate
J = 0.003; %Water Use Efficiency
V = 63; %Advection Rate
R = 100; %Transpiration Rate

Beq = @(P) 0.5*(P*J/M + sqrt((P*J/M)^2-4*L/R)); % Formula for stable, nontrivial 
                                                % equilibrium value of B as a func of P 
kres = 0.001;
kmin = 0;
kmax = 0.5;
kindex = 0;
Pres = 0.1;

sig = zeros(length(241:Pres:300),length(kmin:kres:kmax));
A = zeros(2,2);


for k = kmin:kres:kmax %Wave numbers
    kindex = kindex + 1;
    Pindex = 0;
    for P = 241:Pres:300 %Turing Bif occurs in this range according to Gandhi et al plot
        Pindex = Pindex + 1;
        B = Beq(P);
        W = M/(J*R*B);
        A(1,1) = V*1i*k - R*B^2 - L;
        A(1,2) = -2*R*W*B;
        A(2,1) = J*R*B^2;
        A(2,2) = -D*k^2 + 2*J*R*W*B - M;
        sig(Pindex,kindex) = max(real(eig(A)));
%         if sig(Pindex,kindex) < 0
%            sig(Pindex,kindex) = 0;
%         else
%             sig(Pindex,kindex) = 1;
%         end
    end
end

%plot

xax = kmin:kres:kmax;
yax = 241:Pres:300;

[X, Y] = meshgrid(xax,yax);

contour(X,Y,sig,[0 0])
title('Growing/Decaying Modes of Homogenous Equilibrium Perturbations')
xlabel('Wavenumber')
ylabel('P')
legend('\sigma_{max} = 0')


