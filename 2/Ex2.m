%Import data

dt0 = -7.661711424589D-05;
dt1 = -3.183231456205D-12;
dt2 =  0.000000000000D+00;
sqrta = 5.153650835037D+03; % (sqrt of meters)
a = sqrta^2;
e = 3.841053112410D-03;
M0 = 1.295004883409D+00; %(radians)
Omega0 = -2.241692424630D-01; %(radians)
Omegadot = -8.386063598924D-09; %(radians/sec)
i0 = 9.634782624741D-01; %(radians)
idot = -7.286017777600D-11; %(radians/sec)
w0 = 9.419793734505D-01; %(radians)
wdot = 0.0; %(radians/sec)
GMe = 3.986005D+14; %(m3/s2)
OmegaEdot = 7.2921151467D-05; %(radians)

%CLOCK OFFSET
t0 = 0;
tend = 86399;
tstep = 30;

counter = 0;

for Dt = t0 : tstep : (tend - t0)
    counter = counter+1;
    clockOffsets(1,counter) = Dt;
    clockOffsets(2,counter) = dt0 + dt1*Dt + dt2*Dt^2;
end

% Plot the data
timePlot = plot(clockOffsets(1, :), clockOffsets(2, :));
xlabel('Time (s)');
ylabel('Offset');
title('Clock Offset vs. Time');
ax = ancestor(timePlot, 'axes');
ax.XAxis.Exponent = 0;

fileID = fopen('matlab_results.txt', 'w');

ORSCoord = zeros (3, counter);
ITRFCoord = zeros(3, counter);
ITRFGeo = zeros(3, counter);
epoch = 0;
for Dt = t0 : tstep : (tend - t0)
    epoch = epoch+1;
    n = sqrt(GMe/(a^3));
    Mt = M0 + n*(Dt-t0);
    % compute psi
    eta = estimateEta(0, 10D-10, Mt, e);
    psi = atan2((sqrt(1-e^2)*sin(eta)),(cos(eta)-e));
    % compute r 
    r = (a*(1-e^2))/(1+e*cos(psi));
    % compute w(t)=w0+wdot(t-t0)
    w = w0 + wdot*(Dt-t0);
    % compute i(t)
    i = i0 + idot*(Dt-t0);
    % compute W(t)
    Omega = Omega0 + (Omegadot - OmegaEdot) * (Dt - t0);
    
    % compute x and y
    x = r * cos(psi);
    y = r * sin(psi);
    z = 0;
    ORSCoord(1, epoch) = x;
    ORSCoord(2, epoch) = y;
    ORSCoord(3, epoch) = z;

    %rotate from OCRS to ITRF
    
    ITRFCoord(:,epoch) = ORS2ITRF(ORSCoord(:,epoch), Omega, i, w);   
    ITRFGeo(:,epoch) = gc2gg(ITRFCoord(:,epoch));
    for i = 1:2
    ITRFGeo(i, epoch) = rad2dec(ITRFGeo(i, epoch));
    end

    %Print stuff on the file
    fprintf(fileID, 'Coordinates of the satellite at time: %d\n', Dt);
    fprintf(fileID, 'ORS: X: %.3f Y: %.3f\n', ORSCoord(1, epoch), ORSCoord(2, epoch));
    fprintf(fileID, 'ITRF Cartesian: X: %.3f Y: %.3f Z: %.3f\n', ITRFCoord(1,epoch),  ITRFCoord(2,epoch),  ITRFCoord(3,epoch));
    fprintf(fileID, 'ITRF Geodetic: lat: %.3f Lon: %.3f h: %.3f\n\n', ITRFGeo(1,epoch), ITRFGeo(2,epoch), ITRFGeo(3,epoch));

end
fclose(fileID);

% 4) Plot satellite's daily trajectory with basemap
figure(2);

% H = subplot(m,n,p), or subplot(mnp), breaks the Figure window
% into an m-by-n matrix of small axes

% Plot groundtracks
subplot(3,1,1:2);
% axesm Define map axes and set map properties
ax = axesm ('eqdcylin', 'Frame', 'on', 'Grid', 'on', 'LabelUnits', 'degrees', 'MeridianLabel', 'on', 'ParallelLabel', 'on', 'MLabelParallel', 'south');
% geoshow Display map latitude and longitude data 
%  DISPLAYTYPE can be 'point', 'line', or 'polygon' and defaults to 'line'
geoshow('landareas.shp', 'FaceColor', 'black');
hold on
geoshow(ITRFGeo(1,:),ITRFGeo(2,:), 'DisplayType', 'point', 'MarkerEdgeColor', 'green');
title('Satellite daily trajectory with basemap');
% axis EQUAL  sets the aspect ratio so that equal tick mark
% increments on the x-,y- and z-axis are equal in size.
% axis TIGHT  sets the axis limits to the range of the data.
axis equal; axis tight;

% Plot height of the satellite 
subplot(3,1,3);
plot(clockOffsets(1, :), (ITRFGeo(3,:)-mean(ITRFGeo(3,:)))*0.001, '.g');
title(['Ellipsoidic height variations [km] around mean height = ' num2str(mean(ITRFGeo(3,:))*0.001) ' km']);
xlabel('seconds in one day (00:00 - 23:59 = 86400 sec)');
ylabel('[km]');
xlim([t0, tend]);






