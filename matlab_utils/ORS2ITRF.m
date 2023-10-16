function [itrf] = ORS2ITRF(ors, Omega, i, w)
%ORS2ITRF Summary of this function goes here
%   Detailed explanation goes here

R31 = [cos(-Omega) sin(-Omega) 0; -sin(-Omega) cos(-Omega) 0; 0 0 1];
R1 = [1 0 0; 0 cos(-i) -sin(-i); 0 sin(-i) cos(-i)];
R32 = [cos(-w) sin(-w) 0; -sin(-w) cos(-w) 0; 0 0 1];

itrf = R31*R1*R32*ors;



end

