%Up until the comment line is an excerpt of runAndPlot.m

%Note, this code takes about 3 minutes to run when Pax is 272:-1:100

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

%Initial perturbation is a sine wave with 'numwave' waves in the domain
numwave = 4;
for j = 1:xmax
    Btmax(j) = Beq + 0.1*(sin(2*pi*numwave*j/xmax)); % Initial vegetation density
    Wtmax(j) = Weq + 0.1*(sin(2*pi*numwave*j/xmax)); % Initial water concentration
end

%Preallocations for naive stepping through Turing Equilibria
Pindex = 0;
Pax = 272:-1:150; 
Bavg = zeros(1,length(Pax));
Bmax = Bavg;
numPeaks = Bavg;
wavePeriod = Bavg;
domainPeriod = Bavg;
waveSpeed = Bavg;
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
    BtmaxToPlot(:,Pindex) = Btmax;
    

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
    %Find the location of a peak
    if (numPeaks(Pindex) > 0)
        for i = 2:xmax-1
            if Btmax(i) > Btmax(i-1) && Btmax(i) > Btmax(i+1)
                peakIndex = i;
                BtmaxToPlot(:,Pindex) = cycle(BtmaxToPlot(:,Pindex),peakIndex);
                break
            end
        end
    end

    
    %Now track specifically the point B(peakIndex) and see how many
    %timesteps it will take to reach the peak value again to find the
    %period, and calculate the speed from that.
    BwaveCrest = Btmax(peakIndex);

    %Need the following for KBMEMat_Init to run (These are defaults in f_oneDPDE.m)
    dt = 0.01; %Timestep size
    dx = 0.1; %Grid resolution

    %Backward Euler Matrices
    [lambdaMatInv, muMatInv] = KMBEMat_Init(dt,dx,V,D,xmax);
    
    %Initial timestep and period (in years) timing
    [W_old, B_old] = klausmeierBackwardEuler(Wtmax,Btmax,L,R,J,P,M,dt,lambdaMatInv,muMatInv);
    [W_new, B_new] = klausmeierBackwardEuler(W_old,B_old,L,R,J,P,M,dt,lambdaMatInv,muMatInv);
    W_diff_old = W_old(peakIndex) - W_new(peakIndex);
    W_diff_new = W_diff_old;
    Period = dt;


    %Subsequent timesteps, stop when B(peakIndex) space reaches a peak
    %again 

    %Search for trough
    while W_diff_new*W_diff_old > 0 %while signs are the same
          W_diff_old = W_diff_new;
          W_old = W_new;
          B_old = B_new;
          [W_new, B_new] = klausmeierBackwardEuler(W_old,B_old,L,R,J,P,M,dt,lambdaMatInv,muMatInv);
          W_diff_new = W_old(peakIndex) - W_new(peakIndex);
          Period = Period + dt;
    end

    [W_new, B_new] = klausmeierBackwardEuler(W_old,B_old,L,R,J,P,M,dt,lambdaMatInv,muMatInv);
    W_diff_old = W_old(peakIndex) - W_new(peakIndex);
    W_diff_new = W_diff_old;
    Period = Period + dt;

    while W_diff_new*W_diff_old > 0 %Then check for a new peak
          W_diff_old = W_diff_new;
          W_old = W_new;
          B_old = B_new;
          [W_new, B_new] = klausmeierBackwardEuler(W_old,B_old,L,R,J,P,M,dt,lambdaMatInv,muMatInv);
          W_diff_new = W_old(peakIndex) - W_new(peakIndex);
          Period = Period + dt;
    end
    wavePeriod(Pindex) = Period; %The period of each wavelength calculated by the above algorithm
    domainPeriod(Pindex) = Period*numPeaks(Pindex); %The period of the whole domain: this seems to be what matcont calculates.
    waveSpeed(Pindex) = domlength/domainPeriod(Pindex); %(Phase velocity)
end


%Plot uniform and Turing equilibira
figure(1)
clf
hold on

plot(stable(3,:),stable(2,:), '-k','LineWidth',3);
plot(unstable(3,:),unstable(2,:), '--k','LineWidth',3)
title('Bifurcation Curve for Uniform Vegetation Patterns','FontSize',24)
ylabel('Biomass Density [kgm^{-2}]','FontSize',24)
xlabel('Mean Annual Precipitation [mmyr^{-1}]','FontSize',24)
ax=gca;
ax.FontSize = 24;

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

%Another type of figure
figure(2)
clf
[X,Y] = meshgrid(xax,Pvec);
surf(X,Y,(BtmaxToPlot)','EdgeColor','none')
colormap(gca,flipud(summer))
c = colorbar;
c.Label.String = 'Biomass Density [kgm^{-2}]';
xlabel('Spatial Position [m]')
ylabel('Annual Precipitation [mmyr^{-1}]')
title('Vegetation Density in x-P space')

%Another type of figure
figure(3)
clf
plot(Pvec,waveSpeed)
ylabel('Wave Speed [myr\^{-1}]')
xlabel('Annual Precipitation [mmyr\^{-1}]')
title('Wave speeds')