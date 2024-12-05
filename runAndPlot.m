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


%Plot Precipiation against Biomass with stabilities
figure(1)
clf
hold on
plot(stable(3,:),stable(2,:), '-k','LineWidth',3);
plot(unstable(3,:),unstable(2,:), '--k','LineWidth',3)
xlim([220 320]);
legend('stable','unstable',location = 'northwest')
title('Bifurcation Curve for Uniform Vegetation Patterns','FontSize',24)
ylabel('Biomass Density [kgm^{-2}]','FontSize',24)
xlabel('Mean Annual Precipitation [mmyr^{-1}]','FontSize',24)
ax=gca;
ax.FontSize = 24;
