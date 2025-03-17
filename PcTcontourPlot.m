%Load and plot the P-c-T contours. Make sure vegPat245 is the wd.1
figure(1)
clf
hold on
PcTMat = []; %Initialise PcT Matrix as an empty matrix
Tindex = 0;

for T = 16:4:60
    Tindex = Tindex + 1;
    load(fullfile('matContData',strcat('P-c-T=',num2str(T),'.mat'))) %load data

    [~,minPIndex] = min(params(2,:)); %Find minimum point for plotting curve
    periodFold(:,Tindex+1) = [params(2,minPIndex); params(1,minPIndex)];

    plot(params(2,:),params(1,:))
end
load(fullfile('matContData','fullHopfCurve.mat'))
precipHopfLower = precipHopf(1:1500);
[~,minPIndex] = min(precipHopfLower); %Find minimum point for plotting curve
periodFold(:,1) = [precipHopf(minPIndex); wavespeedHopf(minPIndex)];
xq = linspace(periodFold(1,1),periodFold(1,Tindex+1),100); %Values to evaluate interpolation at
%Spline interpolation for curve of folds of PO contours
periodFoldCurve = interp1(periodFold(1,:),periodFold(2,:),xq,'spline');
plot(precipHopf,wavespeedHopf,'ok')
plot(xq,periodFoldCurve,'*')
load('otherData\420StepDownWaveSpeeds.mat');
plot(Pvec,waveSpeed,'square')
ylabel('Wave Speed [myr\^{-1}]')
xlabel('Annual Precipitation [mmyr\^{-1}]')
title('Periodic Orbits')
legend('Period Contours','','','','','','','','','','','','Hopf Curve','Curve of Folds','Example Sims')
xlim([150 300])
ylim([0 10])