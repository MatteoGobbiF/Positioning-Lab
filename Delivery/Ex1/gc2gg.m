function [gg] = gc2gg(gc, a, b, e)

if nargin == 1 %If a and b and e are not specified, use GRS80

    a = 6378137;
    b = 6356752.3141; 
    f = 1/298.257222100882711243;
    e = sqrt(f*(2-f));
    
    elseif nargin < 4
            
            error('computeRn requires either 1 or 4 input arguments');
end

e2b = (a^2-b^2)/b^2;
r = sqrt(gc(1)^2+gc(2)^2);
psi = atan2(gc(3),(r*sqrt(1-e^2)));

gg(1) = atan2(gc(3)+e2b*b*sin(psi)^3, r-e^2*a*cos(psi)^3);
gg(2) = atan2(gc(2),gc(1));

Rn = computeRn(gg(1), a, e);

gg(3) = r/cos(gg(1))-Rn;




end

