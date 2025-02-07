%Returns a vector of nPoints evenly spaced points between input1 and input2 exclusive
function [output] = evenSpacedInsert(input1,input2,nPoints) 

%Initalising output
output = zeros(nPoints+2,1); 
output(nPoints+2) = input2;
output(1) = input1;

spacing = (input2-input1)/(nPoints + 1); %Spacing between each point

%Evenly spaced points between the 2 inputs
for i = 2:nPoints+1
    output(i) = input1 + (i-1)*spacing;
end
