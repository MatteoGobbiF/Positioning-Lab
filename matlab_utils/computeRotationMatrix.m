function [R] = computeRotationMatrix(Xi, Eta, Alpha)
%Compute the final rotation matrix from the three angles

Rx = computeRx(Xi);
Ry = computeRy(Eta);
Rz = computeRz(Alpha);

R = Rz*Ry*Rx;

end

