load("control_points_db.txt");
load("user_db.txt");

differences = zeros(36, 10);
user_coordinates = zeros(25, 2);
for i = 1:length(user_db)
    for j = 1:length(control_points_db)
        differences(j, :) = (control_points_db(j, 4:13) - user_db(i, 2:11)).^2;
    end
    sumsVector = sum(differences, 2);
    [ minValue, cp_id ] = min(sumsVector);
    user_coordinates(i, 1:2)=control_points_db(cp_id, 2:3);    
end

figure
title("User path")
rectangle("Position", [0 0, 20 10]);
axis([-5 25, -5 15]);
hold on
plot(control_points_db(:, 2), control_points_db(:,3), ".black");
hold on
plot(user_coordinates(:,1), user_coordinates(:,2), "m");

