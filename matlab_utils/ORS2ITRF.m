function [itrf] = ORS2ITRF(ors, Omega, i, w)
%ORS2ITRF Summary of this function goes here
%   Detailed explanation goes here

itrf = computeRz(-Omega)*computeRx(-i)*computeRz(-w)*ors;

end

