
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

sig = zeros(2,60);
A = zeros(2,2);

for k = -1:0.4:1 %Wave numbers
    for P = 241:300 %Turing Bif occurs in this range according to Gandhi et al plot
        B = Beq(P);
        W = M/(J*R*B);
        A(1,1) = V*1i*k - R*B^2 - L;
        A(1,2) = -2*R*W*B;
        A(2,1) = J*R*B^2;
        A(2,2) = -D*k^2 + 2*J*R*W*B - M;
        sig(:,P-240) = eig(A);
    end
    plot(241:300,real(sig(2,:)))
end
