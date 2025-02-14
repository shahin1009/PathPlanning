% Define the rhombus obstacle and start/goal points
clear all; close all;

% Define start and goal points
start_point = [-1, 0.2];
goal_point = [1, 0];

% Number of intermediate points
n_points =3;  % This will create n_points intermediate points

% Define rhombus vertices (centered at origin)
rhombus_width = 1;
rhombus_height = 1;
r = 1;  


% Plot the setup
figure;
hold on;

% Plot rhombus
rhombus_x = [rhombus_width/2, 0, -rhombus_width/2, 0, rhombus_width/2];
rhombus_y = [0, rhombus_height/2, 0, -rhombus_height/2, 0];
fill(rhombus_x, rhombus_y, 'r', 'FaceAlpha', 0.3);

% Plot start and goal points
plot(start_point(1), start_point(2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
plot(goal_point(1), goal_point(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
grid on;
axis equal;
title('Path Planning around Rhombus Obstacle');
xlabel('X-axis');
ylabel('Y-axis');

options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'sqp', ...
    'OptimalityTolerance',10^-12 ...
    ,'StepTolerance',10^-12);

x0 = zeros(2 * n_points, 1);
for j = 1:n_points
    alpha = j / (n_points + 1);
    x0(2*j - 1) = start_point(1) + alpha * (goal_point(1) - start_point(1));
    x0(2*j) = -1;  
end

% Run optimization
[x_opt, fval] = fmincon(@(x) objective(x, start_point, goal_point), x0, [], [], [], [], ...
    repmat([-5, -5], n_points, 1), repmat([5, 5], n_points, 1), @(x) nonlinConstraints(x,r), options);

% Create full path
x_init=reshape(x0, 2, [])';
optimal_points = reshape(x_opt, 2, [])';
optimal_path = [start_point; optimal_points; goal_point];

% Plot the optimal path
plot(optimal_path(:,1), optimal_path(:,2), 'b-', 'LineWidth', 2);
plot(optimal_points(:,1), optimal_points(:,2), 'bo', 'MarkerSize', 5);
plot(x_init(:,1), x_init(:,2), 'ro', 'MarkerSize', 5);


% Display results
fprintf('Total path length: %.2f\n', fval);


 

% plot circle ( the line should be outside of the circle)
theta = linspace(0, 2*pi, 100);


x_circle = 0 + r/2 * cos(theta);
y_circle = 0 + r/2 * sin(theta);

plot(x_circle, y_circle, 'm-', 'LineWidth', 2);

legend('Obstacle', 'Start', 'Goal', 'Path', 'Intermediate Points','initial guess','constraint');




function [c, ceq] = nonlinConstraints(x,r)
    points = reshape(x, 2, [])';
    n_points = size(points, 1);
    
    c = zeros(n_points, 1);
    
    for i = 1:n_points
        x_coord = points(i, 1);
        y_coord = points(i, 2);
        
        c(i) = (r/2)^2 - (abs(x_coord).^2 + abs(y_coord).^2);
        % c(i) = (r/2) - (abs(x_coord) + abs(y_coord));
    end
    
    ceq = [];
end


function total_length_sq = objective(x, start_point, goal_point)

    points = reshape(x, 2, [])';
    

    all_points = [start_point; points; goal_point];
    
    differences = diff(all_points);
    total_length_sq = sum(sum(differences.^2));
end



