function [ sigma_von ] = von_mises( shear_stress, bending_stress)
%VON_MISES Summary of this function goes here
%   Detailed explanation goes here

    sigma1 = bending_stress/2 + sqrt(.5*bending_stress.^2 + shear_stress.^2);
    sigma2 = bending_stress/2 - sqrt(.5*bending_stress.^2 + shear_stress.^2);
    
    sigma_von = sqrt((sigma1.^2 + sigma2.^2 + (sigma1 - sigma2).^2)/2);
end

