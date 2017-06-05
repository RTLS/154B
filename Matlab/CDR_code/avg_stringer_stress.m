function [ stringer_bending_stress ] = avg_stringer_stress( bending_stress, ribLocations )
%AVG_STRINGER_STRESS Summary of this function goes here

numRibs = length(ribLocations);
numStringers = length(bending_stress(:,1));
stringer_bending_stress = zeros(numStringers, (numRibs-1));
for i = 1:(numRibs-1)
    for k = 1:numStringers
        stringer_bending_stress(k,i) = mean([bending_stress(k,max(1,floor(ribLocations(i)))) bending_stress(k,floor(ribLocations(i+1)))]);
    end
end


end

