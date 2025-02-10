%This code should run and plot the 4 to 2 wave collapse
load('otherData\4WaveTo2WaveCollapse.mat')

%We initialise at the 4 wave solution, but with P such that is can only
%support a 2 wave solution
[~,~,~,~,Bmat] = f_oneDPDE(precipAfterCollapse,biomassBeforeCollapse,waterBeforeCollapse); %f_oneDPDE runs for 2000 timesteps, which is 20 years

figure(2)
clf
[X,Y] = meshgrid(xax,0.01:0.01:20); %time mesh is hard coded here
surf(X,Y,Bmat,'EdgeColor','none')
colormap(gca,flipud(summer))
c = colorbar;
c.Label.String = 'Biomass Density [kgm^{-2}]';
xlabel('Spatial Position [m]')
ylabel('Time [yr]')
ylim([0 10])
title('Wave Collapse Transient Dynamics')