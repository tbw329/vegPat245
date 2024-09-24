function x=homKlausmeierPvar(y) 
%Params used for diagram in paper (P treated as variable):
L = 4; %yr^-1 - Evaporation Rate
M = 1.8; %yr^-1 Mortality Rate
J = 0.003; %Ratio of dry to wet mass 
R = 100; %kg^-1yr^-1m^-2 Transpiration Rate
%Func:
x = [y(3) - L*y(1) - R*y(1)*y(2)^2;J*R*y(1)*y(2)^2 - M*y(2)];
end
