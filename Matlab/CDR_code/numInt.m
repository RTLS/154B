function [ integral ] = numInt( dist )
%NUMINT Summary of this function goes here
%   Detailed explanation goes here

    integral = zeros(1, length(dist));
    integral(1) = dist(1);
    for i = 2:length(integral)
       integral(i) = integral(i-1) + dist(i);
    end
end

