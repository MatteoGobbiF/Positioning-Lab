function [eta] = estimateEta(prev, threshhold, Mt, e)
%ESTIMATEETA Summary of this function goes here
%   Detailed explanation goes here
newEta = Mt+e*sin(prev);
if(abs(newEta-prev)<threshhold)
    eta = newEta;
else
    eta = estimateEta(newEta, threshhold, Mt, e);
end

