function [Rx] = computeRx(a)
%Rotation matrix around x axis

Rx = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];

end

