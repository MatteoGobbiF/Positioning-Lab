%getting obs values
run('Ex04_variables.m')

%initialization of parameters

%initialization of parameters
%max number of iterations
max_iterations = 10;

%threshold for convergence
threshold = 0.1;

%approx coordinates of the receiver and clock offset
X_init = [0 0 0 0]'; 

% storage of iterate results
x_estimated = zeros(4,max_iterations);

% Least squares iteration
for i=1:max_iterations

    x_estimate = X_init(1:3);
    % topocentric positions of all satellites
    [azimuth, elevation, ro] = topocent(x_estimate, xyz_sat);

    % approximate geodetic coordinates of the receiver
    x_gg = gc2gg(x_estimate); 
    %x_gg = [phi, labda, h]

    % tropospheric and ionospheric corrections
    tropo_eff = tropo_correction(x_gg(3), elevation);
    iono_eff = iono_correction(x_gg(1), x_gg(2), azimuth, elevation, time_rx, ionoparams);

    % LS known term
    b = ro - s_light*dtS + tropo_eff + iono_eff;

    % LS A matrix
    A = [(x_estimate(1)-xyz_sat(:,1))./ro (x_estimate(2)-xyz_sat(:,2))./ro (x_estimate(3)-xyz_sat(:,3))./ro ones(11,1)];
    N = A' * A;

    % Least square solution for the corrections to the apriori
    obs_corr = pr_C1 - b;
    est_corr = inv(N) * A' * obs_corr;

    % Estimated coordinates of the receiver: 
    % approximate + estimated correction
    x_estimated(:,i) = X_init + est_corr;
    X_init = x_estimated(:,i);

    %check convergence of the result and, in case exit
    if max(abs(est_corr(1:3))) < threshold
        break
    end
    % check at the end that convergence did not fail
    if i == max_iterations
        disp('Convergence failed');
    end
end

% final estimated unknowns
final_coordinates = x_estimated(1:3,i);

% LS residuals and sigma2
% covariance matrix of the estimated coordinates
Qxx = inv(N);
Qx = Qxx(1:3,1:3);

% Rotate and PDOP
R = computeR0(x_gg(1), x_gg(2));
Q = R*Qx*R';
PDOP = sqrt(Q(1,1)+ Q(2,2) + Q(3,3));
dtr = (x_estimated(4,i)-x_estimated(4,i-1))/s_light;

% print results
i_print = sprintf('The total number of iterations is: %d', i);
disp(i_print);
coord_print = sprintf('The coordinates of the receiver are: [%d %d %d]', final_coordinates);
disp(coord_print);
offset_print = sprintf('The clock offset of the receiver is: %d', dtr);
disp(offset_print);
PDOP_print = sprintf('PDOP value is %f', PDOP);
disp(PDOP_print);

%%

%% repeat with CutOfAngle

CutOfAngle = 5;
n_sat = length(dtS);
temp = zeros(n_sat,1);

% extract satellites above cut off
for j=1:n_sat
    if elevation(j) > CutOfAngle
        % store coordinates, otherwise exclude
        temp(j)=elevation(j);
    end
end
[index] = find(temp);
newsat = length(index);

% repeat same computations
n_x_init = [final_coordinates; dtr];
n_x_estimated = zeros(4,max_iterations);

for n_i=1:max_iterations
    n_x_estimate = n_x_init(1:3);
    [n_azimuth, n_elevation, n_ro] = topocent(n_x_estimate, xyz_sat(index,:));    
    n_x_gg = gc2gg(n_x_estimate); 
    %n_x_gg = [n_phi, n_labda, n_h]
    n_x_gg(1:2) = n_x_gg(1:2)*180/pi;
    n_tropo_eff = tropo_correction(n_x_gg(3), n_elevation);
    n_iono_eff = iono_correction(n_x_gg(1), n_x_gg(2), n_azimuth, n_elevation, time_rx, ionoparams);
    n_b = n_ro - s_light*dtS(index) + n_tropo_eff + n_iono_eff;
    n_A = [(n_x_estimate(1)-xyz_sat(index,1))./n_ro (n_x_estimate(2)-xyz_sat(index,2))./n_ro (n_x_estimate(3)-xyz_sat(index,3))./n_ro ones(newsat,1)];
    n_N = n_A' * n_A; 
    n_obs_corr = pr_C1(index) - n_b;
    n_est_corr = inv(n_N) * n_A' * n_obs_corr;
    n_x_estimated(:,i) = n_x_init + n_est_corr;
    X_init = x_estimated(:,i);
    if max(abs(n_est_corr(1:3))) < threshold
        break
    end
    if i == max_iterations
        disp('Cut of Angle Convergence failed');
    end
end
n_final_coordinates = n_x_estimated(1:3,i);
n_Qxx = inv(n_N);
n_Qx = n_Qxx(1:3,1:3);
n_R = computeR0(n_x_gg(1), n_x_gg(2));
n_Q = n_R*n_Qx*n_R';
n_PDOP = sqrt(n_Q(1,1)+ n_Q(2,2) + n_Q(3,3));
n_dtr = (n_x_estimated(4,i)-n_x_estimated(4,i-1))/s_light;

% print
disp('-------------------------------------------------------------');
COA_print = sprintf('Cut off Angle: %d',CutOfAngle);
disp(COA_print);
i_print = sprintf('The total number of iterations (COA) is: %d', i);
disp(i_print);
coord_print = sprintf('The coordinates of the receiver (COA) are: [%d %d %d]', n_final_coordinates);
disp(coord_print);
offset_print = sprintf('The clock offset of the receiver (COA) is: %d', n_dtr);
disp(offset_print);
PDOP_print = sprintf('PDOP (COA) value is %f', n_PDOP);
disp(PDOP_print);
