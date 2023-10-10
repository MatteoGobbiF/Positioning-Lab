function [R0] = computeR0(phi, lambda)
%COMPUTE R0 

R0 = [-sin(lambda) cos(lambda) 0;
    -sin(phi)*cos(lambda) -sin(phi)*sin(lambda) cos(phi);
    cos(phi)*cos(lambda) cos(phi)*cos(lambda) sin(phi)];
end

