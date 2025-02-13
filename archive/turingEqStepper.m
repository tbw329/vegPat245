%Up until the comment line is an excerpt of runAndPlot.m

%Note, this code takes about 3 minutes to run when Pstep is -1

%Params for curve tracking alg
stepsize = 0.2;
nmax = 1000;
%These params are for rootfinding
tol = 1e-5;
maxit = 1000;
%This is for the Jacobian
jach = 1e-5;
%Run curve tracking alg, initial guess for eq and tangent are such that it
%emulates paper best
WBPList = mSciTrackCurve(@homKlausmeierPvar,[0.15;0.15;320],[0;0;-1],stepsize,nmax,tol,maxit,jach);

%Get stabilities of each equilibrium
%Initialise indices for stab/unstab vecs
stabno = 1;
unstabno = 1;

for i = 1:nmax
    J = mSciJacobian(@homKlausmeierPvar, WBPList(:,i),jach);
    J = J(1:2,1:2); %Ignore parts of jacobian dep on P
    eigvals = eig(J);
    if all(real(eigvals) < 0)
        stable(:,stabno) = WBPList(:,i);
        stabno = stabno + 1;
    else
        unstable(:,unstabno) = WBPList(:,i);
        unstabno = unstabno + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Start with P = 272
P = 272;

%Other Params (As in Gandhi et al)
D = 10; %Seed dispersal rate - Something funky going on here
L = 4; %Evap Rate
M = 1.8; %Mortality Rate
J = 0.003; %Water Use Efficiency
V = 63; %Advection Rate
R = 100; %Transpiration Rate

%Initial IC

xmax = 1000;
domlength = 100;

Beq = 0.5*(P*J/M + sqrt((P*J/M)^2-4*L/R)); 
Weq = M/(J*R*Beq);

Wtmax=zeros(xmax,1);
Btmax=Wtmax;

for j = 1:xmax
    Btmax(j) = Beq + 0.1*(sin(0.25*j*domlength/xmax)); % Initial vegetation density
    Wtmax(j) = Weq + 0.1*(sin(0.25*j*domlength/xmax)); % Initial water concentration
end

Pindex = 0;
Pax = 272:-1:100; 
Bavg = zeros(1,length(Pax));
Bmax = Bavg;
numPeaks = Bavg;
for P = Pax
    Pindex = Pindex + 1;
    P
    [Bavg(Pindex),Btmax,Wtmax,xax] = f_oneDPDE(P,Btmax,Wtmax);
    Bmax(Pindex) = max(Btmax);
    %This next part finds the number of peaks in Btmax.
    %Finds difference between each point, finds where they change sign by
    %multiplying adjacent differences.
    diffs = zeros(xmax,1);
    changes = diffs;
    diffs(1) = Btmax(1) - Btmax(xmax);
    for i = 2:xmax
        diffs(i) = Btmax(i) - Btmax(i-1);
        changes(i-1) = diffs(i)*diffs(i-1);
    end
    changes(xmax) = diffs(1)*diffs(xmax);
    numPeaks(Pindex) = (sum(changes<0)/2); %Divide by two so we aren't countring troughs as well. 
end

%Plot uniform and Turing equilibira
figure(1)
clf
hold on
plot(stable(3,:),stable(2,:), '-k');
plot(unstable(3,:),unstable(2,:), '--k')
scatter(Pax,Bmax,[],numPeaks)
scatter(Pax,Bavg,[],numPeaks,'filled')
colormap(gca,'jet')
c = colorbar;
c.Label.String = 'Number of Wavelengths in Domain';
legend('Homogenous Stable','Homogenous Unstable','Maximum Biomass in Patterns','Average Biomass in Patterns',location = 'northwest')
title('Equilibria of the Klausemeier System')
ylabel('Biomass Density kgm^{-2}')
xlabel('Precipitation mmyr^{-1}')