%Load and plot the P-c-T contours. Make sure vegPat245 is the wd.1
figure(1)
hold on
for T = 16:4:60
    load(fullfile('matContData',strcat('P-c-T=',num2str(T),'.mat')))
    plot(params(2,:),params(1,:))
end
load(fullfile('matContData','HopfCurve.mat'))
plot(precipHopf,wavespeedHopf,'ok')
ylabel('Wave Speed[myr\^{-1}]')
xlabel('Precipitation [mmyr\{^-1}]')