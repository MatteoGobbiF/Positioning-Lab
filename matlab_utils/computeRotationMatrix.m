function [R] = computeRotationMatrix(Xi, Eta, Alpha)
%Compute the final rotation matrix from the three angles

Rx = [1 0 0; 0 cos(Xi) -sin(Xi); 0 sin(Xi) cos(Xi)];
Ry = [cos(Eta) 0 -sin(Eta); 0 1 0; sin(Eta) 0 cos(Eta)];
Rz = [cos(Alpha) sin(Alpha) 0; -sin(Alpha) cos(Alpha) 0; 0 0 1];

R = Rz*Ry*Rx;

end

