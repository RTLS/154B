function [ stringerNames ] = get_stringer_names( numNoseTopStringers, numTopStringers, numBottomStringers, numStringers )
%GET_STRINGER_NAMES Summary of this function goes here
%   Detailed explanation goes here

stringerNames = cell(numStringers, 1);
for i = 1:numStringers
   if i <= numNoseTopStringers
       stringerNames(i) = {['topNoseStringer' int2str(i)]};
   elseif i == numNoseTopStringers + 1
       stringerNames(i) = {['frontSparTopCap']};
   elseif i <= numTopStringers + numNoseTopStringers + 1
       stringerNames(i) = {['topStringer' int2str(i - numNoseTopStringers - 1)]};
   elseif i == numTopStringers + numNoseTopStringers + 2
       stringerNames(i) = {['backSparTopCap']};
   elseif i == numTopStringers + numNoseTopStringers + 3
       stringerNames(i) = {['backSparBottomCap']};
   elseif i <= numTopStringers + numNoseTopStringers + numBottomStringers + 3
       stringerNames(i) = {['bottomStringer' int2str(i - numNoseTopStringers - numTopStringers - 3)]};
   elseif i == numTopStringers + numNoseTopStringers + numBottomStringers + 4
       stringerNames(i) = {['frontSparBottomCap']};
   else
       stringerNames(i) = {['bottomNoseStringer' int2str(i - numNoseTopStringers - numTopStringers - numBottomStringers - 4)]};
   end
end

end

