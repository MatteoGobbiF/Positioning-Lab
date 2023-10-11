%Declare constants
a = 6378137;
f = 1/298.257222100882711243;
e = sqrt(f*(2-f));

%Compute origin global geodetic and global cartesian coordinates
Xogeo = [dec2rad(sex2dec(44,23,24.000)); dec2rad(sex2dec(8,56,20.000)); 70];
Xocart = gg2gc(Xogeo, a, e);

%Insert starting data
Xall = [0;30;0];
Xbll = [0;-30;0];
Xcll = [200;0;0];

Xi = dec2rad(sex2dec(0,0,10.23));
Eta = dec2rad(sex2dec(0,0,9.5));
Alpha = dec2rad(sex2dec(30, 27, 18));

%Compute rotation matrix from local level to local cartesian
Rlc2ll = computeRotationMatrix(Xi, Eta, Alpha);
Rll2lc = Rlc2ll';

%Local Cartesian
Xalc = Rll2lc*Xall
Xblc = Rll2lc*Xbll
Xclc = Rll2lc*Xcll

%Global Cartesian
Xagc = lc2gc(Xalc, Xocart, Xogeo(1), Xogeo(2))
Xbgc = lc2gc(Xblc, Xocart, Xogeo(1), Xogeo(2))
Xcgc = lc2gc(Xclc, Xocart, Xogeo(1), Xogeo(2))

%Conversion to ETRF on the website

Xa_ETRF = [4509854.8133; 709344.7329; 4439228.7611];
Xb_ETRF = [4509885.8305; 709380.3979; 4439191.8018];
Xc_ETRF = [4509773.4717; 709521.8492; 4439282.7125];

%Conversion from ETRF Global Cartesian to Global Geodetic

Xaggrad = gc2gg(Xa_ETRF);
Xbggrad = gc2gg(Xb_ETRF);
Xcggrad = gc2gg(Xc_ETRF);

%Convert to sexagesimal

Xagg = [dec2sex(rad2dec(Xaggrad(1)))'; dec2sex(rad2dec(Xaggrad(2)))']
Xaggh = Xaggrad(3)

Xbgg = [dec2sex(rad2dec(Xbggrad(1)))'; dec2sex(rad2dec(Xbggrad(2)))']
Xbggh = Xbggrad(3)

Xcgg = [dec2sex(rad2dec(Xcggrad(1)))'; dec2sex(rad2dec(Xcggrad(2)))']
Xcggh = Xcggrad(3)

% Covariance Propagation

sigmaA = 0.02;
sigmaB = 0.02;
sigmaC = 0.1;
sigmaO = 0.1;

CllA = diag([sigmaA^2, sigmaA^2, sigmaA^2]);
CllB = diag([sigmaB^2, sigmaB^2, sigmaB^2]);
CllC = diag([sigmaC^2, sigmaC^2, sigmaC^2]);
CggO = diag([sigmaO^2, sigmaO^2, sigmaO^2]);

ClcA = Rll2lc*CllA*Rll2lc';
ClcB = Rll2lc*CllB*Rll2lc';
ClcC = Rll2lc*CllC*Rll2lc';

R0 = computeR0(Xogeo(1), Xogeo(2));

CgcA = CggO + R0'*ClcA*R0;
CgcB = CggO + R0'*ClcB*R0;
CgcC = CggO + R0'*ClcC*R0;

Ca = R0*CgcA*R0';
Cb = R0*CgcB*R0';
Cc = R0*CgcC*R0';

%Get the standard deviation in Lat, Lon, Up in cm for each point

deviationA = 100*(sqrt(diag(Ca)))
deviationB = 100*(sqrt(diag(Cb)))
deviationC = 100*(sqrt(diag(Cc)))

%Print the data on a file

fileID = fopen('matlab_results.txt', 'w');

fprintf(fileID, 'Local cartesian coordinates of points A, B and C in meters:\n');
fprintf(fileID, '\nPoint A:\nE: %.3f\nN: %.3f\nU: %.3f\n', Xalc(1), Xalc(2), Xalc(3));
fprintf(fileID, '\nPoint B:\nE: %.3f\nN: %.3f\nU: %.3f\n', Xblc(1), Xblc(2), Xblc(3));
fprintf(fileID, '\nPoint C:\nE: %.3f\nN: %.3f\nU: %.3f\n', Xclc(1), Xclc(2), Xclc(3));

fprintf(fileID, '\nITRF global cartesian coordinates of points A, B, C\n');
fprintf(fileID, '\nPoint A:\nX: %.3f\nY: %.3f\nZ: %.3f\n', Xagc(1), Xagc(2), Xagc(3));
fprintf(fileID, '\nPoint B:\nX: %.3f\nY: %.3f\nZ: %.3f\n', Xbgc(1), Xbgc(2), Xbgc(3));
fprintf(fileID, '\nPoint C:\nX: %.3f\nY: %.3f\nZ: %.3f\n', Xcgc(1), Xcgc(2), Xcgc(3));

fprintf(fileID, '\nETRF global geodetic coordinates of points A, B, C\n');
fprintf(fileID, '\nPoint A:\nLat: %.3f %.3f %.3f\nLong: %.3f %.3f %.3f\nh: %.3f\n', Xagg(1,1), Xagg(1,2), Xagg(1,3), Xagg(2,1), Xagg(2,2), Xagg(2,3), Xaggh);
fprintf(fileID, '\nPoint B:\nLat: %.3f %.3f %.3f\nLong: %.3f %.3f %.3f\nh: %.3f\n', Xbgg(1,1), Xbgg(1,2), Xbgg(1,3), Xbgg(2,1), Xbgg(2,2), Xbgg(2,3), Xbggh);
fprintf(fileID, '\nPoint C:\nLat: %.3f %.3f %.3f\nLong: %.3f %.3f %.3f\nh: %.3f\n', Xcgg(1,1), Xcgg(1,2), Xcgg(1,3), Xcgg(2,1), Xcgg(2,2), Xcgg(2,3), Xcggh);

fprintf(fileID, '\nStandard deviations of points A, B, C in East, North, Up in cm\n');
fprintf(fileID, '\nPoint A:\nEast: %.1f\nNorth: %.1f\nUp: %.1f\n', deviationA(1), deviationA(2), deviationA(3));
fprintf(fileID, '\nPoint B:\nEast: %.1f\nNorth: %.1f\nUp: %.1f\n', deviationB(1), deviationB(2), deviationB(3));
fprintf(fileID, '\nPoint C:\nEast: %.1f\nNorth: %.1f\nUp: %.1f\n', deviationC(1), deviationC(2), deviationC(3));

fclose(fileID);


