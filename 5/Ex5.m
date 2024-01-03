% Import data in Matlab and graph observed DDs
[newdata] = importdata('CycleSlipsDataSun.txt', '	', 1);
idata = newdata.data;

threshold = 3.8*1e-2; % [m]
lambda = 19*1e-2; % [m]

% For each epoch compute differences between observed DDs and approximated DDs (residual DDs) and graph them

difference = zeros(100,2);
for i=1:100
    difference(i,1) = i;
    difference(i,2)=idata(i,2)-idata(i,3);
end


figure
plot(difference(:,1),difference(:,2))

% Compute differences between consecutive epochs of residual DDs (hint: diff or for cycle) and graph them
difference_epoch(:,2) = diff(difference(:,2));
for i=1:99
    difference_epoch(i,1)=i;
end

figure
plot(difference_epoch(:,1),difference_epoch(:,2))


% Identify cycle slips and repair them (hint: just one for cycle with both the actions)

new_DD = idata;
new_difference = difference;
new_difference_epoch = difference_epoch;
while(true)
    [new_DD, cycle_slip] = correct_cycle_slip(new_DD,new_difference_epoch,lambda,threshold);

    for i=1:100
    new_difference(i,1) = i;
    new_difference(i,2)=new_DD(i,2)-idata(i,3);
    end

    new_difference_epoch(:,2) = diff(new_difference(:,2));
    for i=1:99
        new_difference_epoch(i,1)=i;
    end

    if cycle_slip==false
        break
    end
   
end


figure
plot(new_difference(:,1),new_difference(:,2))

