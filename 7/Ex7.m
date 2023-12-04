clc
clear all

%load data without errors
data = load("Inertial_data.dat");
epoch = data(:, 1);
%accelerations
acc = data(:, 2:3); 
%omega in zed
omegaz = data(:,4);

%load data with errors (repeat)
data_err = load("Inertial_data_ni.dat");
epoch_err = data_err(:, 1);
%accelerations
acc_err = data_err(:, 2:3); 
%omega in zed
omegaz_err = data_err(:,4);

%call function to compute trajectory without errors
traject = CalcTrajectory(acc, omegaz, epoch);
%call function to compute trajectory with errors
traject_err = CalcTrajectory(acc_err, omegaz_err, epoch);
%plot comparison on the same plot

figure
plot(traject(:,1),traject(:,2),'-c');
hold on
plot(traject_err(:,1),traject_err(:,2),'-m');
xlabel('x[m]');
ylabel('y[m]');
legend('Without errors','With errors','Location','southwestoutside');