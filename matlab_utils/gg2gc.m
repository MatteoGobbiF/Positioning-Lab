function [gc] = gg2gc(gg, a, e)
%GG2GC Global geodetic to global cartesian conversion

Rn = computeRn(gg(1), a, e);

x = (Rn + gg(3))*cos(gg(1))*cos(gg(2));
y = (Rn + gg(3))*cos(gg(1))*sin(gg(2));
z = (Rn*(1-e^2) + gg(3))*sin(gg(1)); 

gc = [x;y;z];

end

