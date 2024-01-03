function [dec] = sex2dec(degrees, minutes, seconds)
%Sexagesimal to decimal degrees
dec = degrees + (minutes/60) + seconds/3600;
end

