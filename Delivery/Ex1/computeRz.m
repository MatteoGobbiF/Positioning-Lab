function [Rz] = computeRz(a)
%Rotation matrix around z axis

Rz = [cos(a) sin(a) 0; -sin(a) cos(a) 0; 0 0 1];

end

