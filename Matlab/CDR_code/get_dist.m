function [ distance ] = get_dist( xstart, xend, suface)
%GET_DIST Summary of this function goes here
%   Detailed explanation goes here

    % Right now straight line analysis... FIXME
    zstart = get_z(xstart, surface);
    zend = get_z(xend, surface);
    
    distance = sqrt((zend-zstart)^2 + (xend-xstart)^2);

end

