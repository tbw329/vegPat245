function [vecOut] = cycle(vecIn,n) %Cycles the elements in a column vec st element n moves to position 1
vecOut = zeros(length(vecIn),1); %preallocating vecOut

%This cycles the elements
repVec = [vecIn;vecIn];
index = 0;
for i = n:length(vecIn)+n-1
    index = index + 1;
    vecOut(index) = repVec(i);
end
