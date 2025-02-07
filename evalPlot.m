%Reset Plots
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


%Values of k, P to loop over
kres = 0.001;
kmin = 0;
kmax = 0.5;
kvec = kmin:kres:kmax;
kindex = 0;
Pmin = 241; %Turing Bif occurs in 241<P<300 according to Gandhi et al plot
Pmax = 300;
Pres = 0.1;
Pvec = Pmin:Pres:Pmax;

%Pre-allocating sigma(P,k) matrix and eigenvalue matrix
sig = zeros(length(Pvec),length(kvec));
A = zeros(2,2);


for k = kvec %Wave numbers
    kindex = kindex + 1;
    Pindex = 0; %Reset Pindex for every value of k
    for P = Pvec %Precip. Levels
        Pindex = Pindex + 1;
        B = Beq(P);
        W = M/(J*R*B);
        %Constructing eigenvalue matrix
        A(1,1) = V*1i*k - R*B^2 - L;
        A(1,2) = -2*R*W*B;
        A(2,1) = J*R*B^2;
        A(2,2) = -D*k^2 + 2*J*R*W*B - M;
        %Finding eigenvalues
        sig(Pindex,kindex) = max(real(eig(A)));
    end
end

%plotting

%Fig 1: sigma(k) contour plot in k-P space
[X, Y] = meshgrid(kvec,Pvec);
figure(1)
clf
hold on
fig1 = gca;
contourf(fig1,X,Y,sig);
contour(fig1,X,Y,sig,[0 0],'--r','LineWidth',5);
colormap(gca,'parula')
c = colorbar;
c.Label.String = 'max(\Re(\sigma(k)))';
c.FontSize = 24;
ax=gca;
ax.FontSize = 24;
title('Growth rates Wave-Like Equilibrium Perturbations','FontSize',24)
xlabel('Wavenumber [m^{-1}]','FontSize',24)
ylabel('Mean Annual Precipitation [mmyr^{-1}]','FontSize',24)



%Fig 3: sig(k) plot for a range of P values
figure(3)
clf
fig3 = gca;
pl3 = plot(fig3,kvec,sig(1:100:501,:));
title('max(\Re(\sigma(k))) for Selected Values of P')
xlabel('k')
ylabel('max(\Re(\sigma(k)))')
legend('P = 240','P = 250', 'P = 260', 'P = 270', 'P = 280', 'P = 290', 'location','south')