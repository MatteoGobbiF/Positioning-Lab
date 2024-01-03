%load parameters
line1 = 'GPSA 7.4506D-09 1.4901D-08 -5.9605D-08 -1.1921D-07 IONOSPHERIC CORR';
line2 = 'GPSB 9.2160D+04 1.3107D+05 -6.5536D+04 -5.2429D+05 IONOSPHERIC CORR';
ionoparams = [cell2mat(textscan(line1, '%*s %f %f %f %f %*s')) ...
cell2mat(textscan(line2, '%*s %f %f %f %f %*s'))];

%Zenithal maps
el = 90;
az = 0;

tstep = 3600*6;
gstep = 0.5;
t0 = 0;
tend = tstep*3;
lat0 = -80;
latend = 80;
lon0 = -180;
lonend = 180;

Dt=t0:tstep:tend;
lats=lat0:gstep:latend;
lons=lon0:gstep:lonend;

matrix = zeros(numel(lons), numel(lats), numel(Dt));

for t = 1:length(Dt)
    current_time = Dt(t);  % Update the current time value
    for lon = 1:length(lons)
        for lat = 1:length(lats)
            matrix(lon, lat, t) = iono_correction(lats(lat), lons(lon), az, el, current_time, ionoparams);
        end
    end   
end


%Plot the zenithal map

[phi_grid,lambda_grid] = meshgrid(lats,lons);
figure(1)
title('Ionospheric Error Maps')
for i = 1:length(Dt)
    subplot(2,2,i);
    geoshow(phi_grid, lambda_grid, matrix(:,:,i), 'DisplayType','texturemap','facealpha',.5)
    hold on
    geoshow('landareas.shp', 'FaceColor', 'none');
    title(['time = ', num2str(Dt(i)/3600),':00']);
    xlabel('longitude [deg]')
    ylabel('latitude [deg]')
    xlim([-180 180]);
    ylim([-80 80]);
    colormap(jet);
end
hp4 = get(subplot(2,2,4),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.028  hp4(2)  0.03  hp4(2)+hp4(3)*2.1]);

%Polar map in Milano
%El=0:90 az=-180:180 
% Milano position in degrees
phi = 45 + 28 / 60 + 38.28 / 60^2; %degrees
lambda = 9 + 10 / 60 + 53.40 / 60^2; %degrees

t0 = 0;
tstep = 3600*12;
tend = tstep;
el0=0;
elstep = 0.5;
elend=90;
az0=-180;
azstep=0.5;
azend=180;

Dt=t0:tstep:tend;
els = el0:elstep:elend;
azs = az0:azstep:azend;
matrixPolar = zeros(numel(els), numel(azs),numel(Dt));

for t = 1:length(Dt)
    current_time = Dt(t);  % Update the current time value
    for el = 1:length(els)
        for az = 1:length(azs)
            matrixPolar(el, az, t) = iono_correction(phi, lambda, azs(az), els(el), current_time, ionoparams);
        end
    end   
end

%Polar plot

[Az, El] = meshgrid(azs, els);
for i = 1 : 2
    figure(i + 1)
    title(['Ionospheric Error Polar Map for Milan Observer time = ', num2str(Dt(i)/3600),':00'])
    axesm('eqaazim', 'MapLatLimit', [0 90]);
    axis off
    framem on
    gridm on

    mlabel on
    plabel on;
    setm(gca,'MLabelParallel',0)
    geoshow(El, Az, matrixPolar(:,:,i), 'DisplayType','texturemap', 'facealpha',.6)
    colormap(jet)

    hcb = colorbar('eastoutside');
    set(get(hcb,'Xlabel'),'String','Legend')
end



