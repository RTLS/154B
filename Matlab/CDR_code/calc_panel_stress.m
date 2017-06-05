function [ panel_stress ] = calc_panel_stress( stringer_stress )
%CALC_PANEL_STRESS Summary of this function goes here
%   Detailed explanation goes here
numPanelRows = length(stringer_stress(:,1)) - 1;
numPanelCols = length(stringer_stress(1,:));

panel_stress = zeros(numPanelRows, numPanelCols);

for i = 1:numPanelRows
   for j = 1:numPanelCols
       panel_stress(i,j) = mean([stringer_stress(i,j) stringer_stress(i+1,j)]);
   end
end

end

