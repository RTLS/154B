function [ bending_stress ] = calc_bending_stress(MX, MZ, Ix, Iz, Ixz, locationsX, locationsZ, centroid )
%CALC_BENDING_STRESS Summary of this function goes here
%   Detailed explanation goes here

numStringers = length(locationsX);
numCrossSections = length(MX);

bending_stress = zeros(numStringers, numCrossSections);

for i = 1:numCrossSections
    frac1 = (Ix.*MZ(i) + Ixz.*MX(i))/(Ix.*Iz - Ixz.^2);
    frac2 = (Iz.*MX(i) + Ixz.*MZ(i))/(Ix.*Iz - Ixz.^2);
    bending_stress(:, i) = (locationsX - centroid.posX).*frac1 - (locationsZ - centroid.posZ).*frac2;
end

end

