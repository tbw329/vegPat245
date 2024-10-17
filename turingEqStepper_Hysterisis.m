%Up until the comment line is an excerpt of runAndPlot.m

%Here, we step down until we're on the two-wave state, then step up again

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

%Initial ICs

%Defining space discretisation
xmax = 1000;
domlength = 100;

%Analytical form of the homegenous equilibria
Beq = 0.5*(P*J/M + sqrt((P*J/M)^2-4*L/R)); 
Weq = M/(J*R*Beq);

%Preallocating Wtmax and Btmax
Wtmax=zeros(xmax,1);
Btmax=Wtmax;

%Initial perturbation is a sine wave with wavenumber 4 domlength^-1 
for j = 1:xmax
    Btmax(j) = Beq + 0.1*(sin(2*pi*4*j/xmax)); % Initial vegetation density
    Wtmax(j) = Weq + 0.1*(sin(2*pi*4*j/xmax)); % Initial water concentration
end

%Preallocations for naive stepping through Turing Equilibria
Pindex = 0;
Pax = 272:-1:100; 
Bavg = zeros(1,length(Pax));
Bmax = Bavg;
numPeaks = Bavg;
Pvec = [];

%Initially, P decreases from P = 272

for P = Pax
    
    %Indexing
    Pindex = Pindex + 1;
    Pvec = [Pvec,P]; 
    P
    
    %Actually solving the system from the previous steady state (and finding Bavg and Bmax to plot)
    [Bavg(Pindex),Btmax,Wtmax,xax] = f_oneDPDE(P,Btmax,Wtmax);
    Bmax(Pindex) = max(Btmax);

    %If the amplitude of waves is negligible, we declare system to be
    %spatially homogenous
    if Bmax(Pindex)-Bavg(Pindex) < 1e-3
        numPeaks(Pindex) = 0;
    else
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
        numPeaks(Pindex) = (sum(changes<0)/2); %Divide by two so we aren't counting troughs as well. 
    end
    %Once we join the 2-wavelength branch, break this loop and begin to
    %increase P to observe hysteresis
    if ~(numPeaks(Pindex) == 4)
        newP = P;
        break
    end
end

%Same as before, but now P is increasing to P = 300
for P = newP:1:300

    %Indexing
    Pindex = Pindex + 1;
    Pvec = [Pvec,P];
    P

    %Actually solving the system from the previous steady state (and finding Bavg and Bmax to plot)
    [Bavg(Pindex),Btmax,Wtmax,xax] = f_oneDPDE(P,Btmax,Wtmax);
    Bmax(Pindex) = max(Btmax);
    
    %If the amplitude of waves is negligible, we declare system to be
    %homogenous
    if Bmax(Pindex)-Bavg(Pindex) < 1e-3
        numPeaks(Pindex) = 0;
    else
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
        numPeaks(Pindex) = (sum(changes<0)/2); %Divide by two so we aren't counting troughs as well. 
    end
end


%Plot uniform and Turing equilibira
figure(1)
clf
hold on

%Plotting the Stable and Unstable Steady States
plot(stable(3,:),stable(2,:), '-k');
plot(unstable(3,:),unstable(2,:), '--k')

%Plotting the average and maximum values of the Turing equilibria, coloured
%according to the number of peaks in the domain
scatter(Pvec,Bmax,[],numPeaks)
scatter(Pvec,Bavg,[],numPeaks,'filled')
colormap(gca,'jet')
c = colorbar;
c.Label.String = 'Number of Wavelengths in Domain';

%Legend, titles etc.
legend('Homogenous Stable','Homogenous Unstable','Maximum Biomass in Patterns','Average Biomass in Patterns',location = 'northwest')
title('Equilibria of the Klausemeier System')
ylabel('Biomass Density kgm^{-2}')
xlabel('Precipitation mmyr^{-1}')