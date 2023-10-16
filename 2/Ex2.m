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
    psi = 1/tan((sqrt(1-e^2)*sin(eta))/(cos(eta)-e));
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
end

