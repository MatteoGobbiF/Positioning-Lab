function [Rn] = computeRn(latitude, a, e)
%Compute the gran normal radius
narginchk(1, 3); %Check that the function has between 1 and 3 parameters

if nargin == 1 %If a and e are not specified, use GRS80

    a = 6378137;
    f = 1/298.257222100882711243;
    e = sqrt(f*(2-f)); 
    
    elseif nargin == 2
            
            error('computeRn requires either 1 or 3 input arguments');        

end

Rn = a/sqrt(1 -  (e^2)*sin(latitude)^2); 

end

