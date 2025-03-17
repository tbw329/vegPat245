%Interpolate the MatContData to find the function k(P,c) and plot it
figure(1)
clf
hold on
PcTMat = []; %Initialise PcT Matrix as an empty matrix
for T = 16:4:60
    load(fullfile('matContData',strcat('P-c-T=',num2str(T),'.mat')))
    PcTMat = [PcTMat,[params(2,:);params(1,:);T*ones(1,length(params(2,:)))]];
end
PckInterp = scatteredInterpolant((PcTMat(1:2,:))',((2*pi)./(PcTMat(2,:).*PcTMat(3,:)))','natural','none'); %Wavenumber interpolated surface, using k = 2pi/cT
[PckX,PckY] = meshgrid(1:0.1:7,100:0.1:300); %Mesh of P and c values
Pck = PckInterp(PckY,PckX); %Get a function for k(P,c) to plot
PckContour = surfz(PckY,PckX,Pck);
PcKContour.EdgeColor = 'none';
%PckContour = contour(PckY,PckX,Pck,(2*pi./(100./(1:6))));
% load('otherData\420StepDownWaveSpeeds.mat');
% plot(Pvec,waveSpeed,'+')
% ylabel('Wave Speed [myr\^{-1}]')
% xlabel('Annual Precipitation [mmyr\^{-1}]')
% title('Periodic Orbits')

