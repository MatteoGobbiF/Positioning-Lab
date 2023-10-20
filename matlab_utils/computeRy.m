function [Ry] = computeRy(a)

%   Rotation matrix around y axis
Ry = [cos(a) 0 -sin(a); 0 1 0; sin(a) 0 cos(a)];

end

