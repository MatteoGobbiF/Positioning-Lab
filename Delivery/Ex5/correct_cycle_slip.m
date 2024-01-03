function [new_DD, cycle_slip] = correct_cycle_slip(old_DD, differences, lambda, threshold)
%CORRECT_CYCLE_SLIP Summary of this function goes here
%   Detailed explanation goes here

new_DD = zeros(length(old_DD),2);
new_DD(1,1) = 1;
new_DD(1,2) = old_DD(1,2);

cycle_slip = false;
for i=1:length(differences)
    new_DD(i+1,1) = i+1;
    if abs(differences(i,2))>threshold && ~cycle_slip 
        x = differences(i,2)/lambda;
        n = round(x);
        if lambda * abs(n-x) <= threshold
            disp 'found a cycle slip'
            cycle_slip = true;
        end
    end
    
    if cycle_slip
        new_DD(i+1,2) = old_DD(i+1,2)-lambda*n;
    else
        new_DD(i+1,2) = old_DD(i+1,2);
    end

end

end

