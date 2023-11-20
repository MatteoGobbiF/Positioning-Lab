%% Antenna data
dataset = readtable('S3-T1.csv');

N = table2array(dataset(:, 'Northing'));
E = table2array(dataset(:, 'Easting'));
zone = repmat('32 T', [length(N),1]);

[lat, long] = utm2deg(E, N, zone);

%% Phone data

dataset_m = load('S3-T1.mat');

lat_m = dataset_m.Position.latitude;
long_m = dataset_m.Position.longitude;

[N_m, E_m, zone_m] = deg2utm(lat_m, long_m);

%% Compute distance

stoner = [E, N];
phone = [E_m, N_m];

dist = zeros(length(N_m), 1);



for i = 1:length(N_m)
    Ei = E_m(i);
    Ni = N_m(i);
    distance = sqrt(sum((E-Ei).^2 + (N-Ni).^2, 2));
    dist(i) = min(distance);
end

%mean(dist)
%% Plot stuff

figure;
plot(lat,long,'-r')
hold on
plot(lat_m,long_m,'-b')
