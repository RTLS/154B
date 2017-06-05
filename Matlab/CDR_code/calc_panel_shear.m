function [ panel_shear ] = calc_panel_shear( shear_stress )
%CALC_PANEL_STRESS Summary of this function goes here
%   Detailed explanation goes here
numPanelRows = length(shear_stress(:,1));
numPanelCols = length(shear_stress(1,:));

panel_shear = zeros(numPanelRows, numPanelCols);

panel_shear(1,:) = mean([shear_stress(1,:) (-1*shear_stress(numPanelRows,:))]);
for i = 2:numPanelRows
   for j = 1:numPanelCols
       panel_shear(i,j) = mean([shear_stress(i-1,j) (-1*shear_stress(i,j))]);
   end
end

end