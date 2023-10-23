function [I] = klobuchar(lat, long, E, A, alpha, beta, tgps)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%1. Calculate the earth-centred angle (elevation E in semicircles).
psi = 0.0137/(E + 0.11)-0.022;

%2. Compute the latitude of the Ionospheric Pierce Point (IPP)
lati = lat + psi*cos(A);
if(lati>0.416)
    lati = 0.416;
end
if(lati<-0.416)
    lati = -0.416;
end

%3. Compute the longitude of the IPP.
longi = long + (psi*sin(A)/cos(lati));

%4. Find the geomagnetic latitude of the IPP.
latm = lati + 0.064*cos(longi-1.617);

%5. Find the local time at the IPP.
t = 43200*longi+tgps;
if(t>=86400)
    t=t-86400;
end
if(t<0)
    t=t+86400;
end

%6. Compute the amplitude of ionospheric delay.
for i = 1:3
    
end



end

