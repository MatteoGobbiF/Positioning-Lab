function [gc] = lc2gc(lc,Xogc, phi, lambda)

%Conversion from local cartesian to global cartesian

R = computeR0(phi, lambda);

gc = Xogc + R'*lc;

end

