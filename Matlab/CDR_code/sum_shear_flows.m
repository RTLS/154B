function [ webs ] = sum_shear_flows( webs, qt_X, qt_Y, qt_Z, qs_X, qs_Z )
%SUM_SHEAR_FLOWS Summary of this function goes here
%   Detailed explanation goes here

num = size(webs);
num = num(2);

for i = 1:num
    webs(i).qtot_X = webs(i).qPrime_X + qt_X + qs_X;
    webs(i).qtot_Y = qt_Y;
    webs(i).qtot_Z = webs(i).qPrime_Z + qt_Z+ qs_Z;
end


end

