function [x] = CalcTrajectory(acc, omegaz, epoch)
%CALCTRAJECTORY Summary of this function goes here
%   Detailed explanation goes here

deltat = 1;
alphaold = 0;
vxold = 0;
vyold = 0;
x = zeros(length(epoch), 2);
x(1, :) = [100, 100];


for t=1:length(epoch)
    
    vx = vxold + acc(t, 1)*deltat;
    deltax = vx*deltat + (acc(t, 1)*deltat^2)/2;

    ayclean = acc(t, 2) - omegaz(t)*vx;
    vy = vyold + ayclean*deltat;
    deltay = vy*deltat + (ayclean*deltat^2)/2;
    
    alpha = alphaold + omegaz(t)*deltat;
    R = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];

    delta = [deltax deltay];
    x(t+1, :)=x(t, :)+delta*R';

    
    vxold = vx;
    vyold = vy;
    alphaold = alpha;

end



end

