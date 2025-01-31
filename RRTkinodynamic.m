clear;close all;

% Algorithm Parameters
max_iter = 30000000;
step_size = 10;
goal_radius = 15;
v_max = 4;     
omega_max = 8; 
min_turning_radius = 20.0;  

profile on;
flag_show_runtime = 1;  % choose 1 to show the runtime plot (slow)


% map1
% img1 = imread("Pictures/map_1_d.png"); 
% start= [8,116,0];
% goal = [274,115,pi/2];

% % map2
% img1 = imread("Pictures/map_2_d.png"); 
% start= [9,79,0];
% goal = [262,281,pi/2];
% 
% % map3
% img1 = imread("Pictures/map_3_d.png"); 
% start= [6,6,0];
% goal = [5,298,pi];

% % map4
% img1 = imread("Pictures/map_4_d.png"); 
% start= [11,74,0];
% goal = [191,181,pi/2];

% % map,png
img1 = imread("Pictures/map.pgm");
start= [65,245,3*pi/2];
goal = [214,148,pi/2];

% maze 3
% img1 = imread("Pictures/maze3.jpg");
% start=[85,596,3*pi/2];
% goal = [1056,3,3*pi/2];


img = imbinarize(im2gray(img1), 0.8);
se = strel('disk', 2); 
img = ~imdilate(~img, se);
map_size = size(img);


figure; % Create a figure handle
set(gcf, 'Position', [100, 100, 800, 600]); 
imshow(img, 'InitialMagnification', 'fit'); 
hold on
title('Select START point, then GOAL point');

plotCircle(goal, goal_radius)



%%%%%%%%%%%%%%%%%%% Main algorithm %%%%%%%%%%%%%

tic;
path = kinodynamic_rrt(img, start, goal, max_iter, step_size, goal_radius,v_max,omega_max,min_turning_radius,flag_show_runtime);
elapsedtime = toc;

%% Bottleneck analysis
p = profile('info');
prefix = 'RRTkinodynamic>';
function_table = p.FunctionTable;
filtered_functions = function_table(contains({function_table.FunctionName}, prefix));
for i = 1:length(filtered_functions)
    fprintf('Function: %s\n', filtered_functions(i).FunctionName);
    fprintf('Total Time: %.4f seconds\n', filtered_functions(i).TotalTime);
    fprintf('Number of Calls: %d\n', filtered_functions(i).NumCalls);
    fprintf('-----------------------------\n');
end

%% Functions
function path = kinodynamic_rrt(img, start, goal, max_iter, step_size, goal_radius, v_max, omega_max, min_turning_radius, flag_show_runtime)
    tree.nodes = start;    % [x, y, theta]
    tree.parent = 0;       % Parent index
    tree.cost = 0;         % Cost from start
    tree.controls = [];    % Store control inputs
    tree.trajectories = cell(1); % Store complete trajectories
    tree.path_costs = 0;   % Add path_costs to store distance-based costs

    % Plot start and goal
    plot(start(1), start(2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(goal(1), goal(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    
    % Initialize collision checking parameters
    robot_radius = 2; % Define robot radius for collision checking
    [img_height, img_width] = size(img);
    
    for iter = 1:max_iter
        % 1. Sample a random point with goal biasing
        if rand() < 0.2
            rand_point = goal;
        else
            rand_point = [randi([1, img_width]), randi([1, img_height]), ...
                         (rand() * 2) * pi];
        end

        if mod(iter,1000)==0
        disp(iter);
        end

        % 2. Find nearest node
        [nearest_idx, nearest_node] = find_nearest_node(tree.nodes, rand_point);

        % 3. Generate controls and select best
        best_state = [];
        best_control = [];
        best_trajectory = [];
        best_path_cost = inf;
        min_dist = inf;
        
        for i = 1:4
            control = sample_control(v_max, omega_max, min_turning_radius);
            [new_state, trajectory] = generate_trajectory_euler(nearest_node, control, step_size);
            
            if is_trajectory_collision_free(img, trajectory, robot_radius)
                dist = compute_distance(new_state, rand_point);
                path_cost = compute_trajectory_cost(trajectory);
                
                if dist < min_dist
                    min_dist = dist;
                    best_state = new_state;
                    best_control = control;
                    best_trajectory = trajectory;
                    best_path_cost = path_cost;
                end
            end
        end

        if isempty(best_state)
            continue;
        end

        % 4. Add new node to tree with updated cost calculation
        current_idx = size(tree.nodes, 1) + 1;
        tree.nodes(current_idx, :) = best_state;
        tree.parent(current_idx) = nearest_idx;
        tree.controls(current_idx, :) = best_control;
        tree.path_costs(current_idx) = tree.path_costs(nearest_idx) + best_path_cost;
        tree.trajectories{current_idx} = best_trajectory;

        % Plot the trajectory
        if flag_show_runtime
            plot(best_trajectory(:,1), best_trajectory(:,2), 'b-', 'LineWidth', 0.5);
            % plot(best_state(1), best_state(2), 'b.', 'MarkerSize', 5);
            drawnow limitrate;
        end
        
        % 5. Check if goal reached
        if norm(best_state(1:2) - goal(1:2)) <= goal_radius 
           % && abs(angle_diff(best_state(3), goal(3))) <= pi/8
            [path, final_cost] = reconstruct_path(tree);
            fprintf('Final path cost: %.2f\n', final_cost);
            return;
        end
    end
    
    path = [];
    warning('Path not found within maximum iterations');
end
%% Generate Trajectory

function [new_state, trajectory] = generate_trajectory_euler(state, control, step_size)
    
    dt = (0.5 + 1.5 * rand());
    % dt = 0.1; % Fixed time step for Euler integration
    num_steps = step_size;
    trajectory = zeros(num_steps + 1, 3);
    trajectory(1, :) = state;
    
    current_state = state;

    
    % Get control inputs
    v = control(1);
    omega = control(2);

    if abs(omega)<0.07
        omega=0;
    end
    
    % Generate trajectory using Forward Euler
    for i = 1:num_steps
        % State derivatives
        x_dot = v * cos(current_state(3));
        y_dot = v * sin(current_state(3));
        theta_dot = omega;
        
        % Euler integration
        current_state(1) = current_state(1) + x_dot * dt;
        current_state(2) = current_state(2) + y_dot * dt;
        current_state(3) = current_state(3) + theta_dot * dt;
        
        % Normalize angle
        current_state(3) = mod(current_state(3) + pi, 2*pi) - pi;
        
        trajectory(i + 1, :) = current_state;
    end
    
    new_state = current_state;
end

function cost = compute_trajectory_cost(trajectory)
    % Initialize cost
    cost = 0;
    
    % Calculate Euclidean distance between consecutive points
    for i = 1:size(trajectory,1)-1
        % Get current and next point
        current_point = trajectory(i, 1:2);
        next_point = trajectory(i+1, 1:2);
        
        % Add Euclidean distance to cost
        segment_dist = norm(next_point - current_point);
        cost = cost + segment_dist;
    end
end

function dist = compute_distance(state1, state2)
    pos_diff = norm(state1(1:2) - state2(1:2));
    ang_diff = abs(angle_diff(state1(3), state2(3)));
    dist = pos_diff + 0.1 * ang_diff;
end

function [nearest_idx, nearest_node] = find_nearest_node(nodes, point)
    distances = zeros(size(nodes, 1), 1);
    for i = 1:size(nodes, 1)
        distances(i) = compute_distance(nodes(i,:), point);
    end
    [~, nearest_idx] = min(distances);
    nearest_node = nodes(nearest_idx, :);
end



%%%%%%%%%%%%%%%%%%%%%%%%%% random sampling control %%%%%%%%%%%%%%%%%%
function control = sample_control(v_max,omega_max,min_turning_radius)
    
    % Sample velocity with bias towards higher speeds for better exploration
    v = abs(randn() * (v_max/2));
    v = min(max(v, 0.5), v_max);  % Ensure minimum velocity
    
    % Sample angular velocity
    omega = randn() * (omega_max/3);
    omega = min(max(omega, -omega_max), omega_max);
    
    % Enforce minimum turning radius constraint
    if abs(omega) > 1e-6  % If turning
        radius = abs(v / omega);
        if radius < min_turning_radius
            % Adjust omega to maintain minimum turning radius
            omega = sign(omega) * abs(v / min_turning_radius);
        end
    end
    
    control = [v, omega];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% colision check %%%%%%%%%%%%%%

function valid = is_trajectory_collision_free(img, trajectory, robot_radius)
    valid = true;
    [img_height, img_width] = size(img);
    
    for i = 1:size(trajectory, 1)
        x = round(trajectory(i, 1));
        y = round(trajectory(i, 2));
        
        % Check if point is within image bounds
        if x < 1 || x > img_width || y < 1 || y > img_height
            valid = false;
            break;
        end
        
        % Check collision within robot radius
        for dx = -robot_radius:robot_radius
            for dy = -robot_radius:robot_radius
                if dx^2 + dy^2 <= robot_radius^2
                    check_x = round(x + dx);
                    check_y = round(y + dy);
                    
                    % Ensure checking point is within bounds
                    if check_x >= 1 && check_x <= img_width && ...
                       check_y >= 1 && check_y <= img_height
                        if ~img(check_y, check_x)  % Obstacle detected
                            valid = false;
                            return;
                        end
                    end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%% Reconstruct the path to goal %%%%%%%%%%%%%%%%%%%%%%%%%
function [path, total_cost] = reconstruct_path(tree)
    idx = size(tree.nodes, 1);
    path = tree.nodes(idx, :);
    full_trajectory = tree.trajectories{idx};
    total_cost = tree.path_costs(idx);
    
    while tree.parent(idx) ~= 0
        prev_idx = tree.parent(idx);
        % Plot final path segments
        trajectory = tree.trajectories{idx};
        plot(trajectory(:,1), trajectory(:,2), 'm-', 'LineWidth', 3);
        path = [tree.nodes(prev_idx, :); path];
        full_trajectory = [tree.trajectories{prev_idx}; full_trajectory];
        idx = prev_idx;
    end
    
    path = full_trajectory;
end


function diff = angle_diff(a1, a2)
    diff = mod(a1 - a2 + pi, 2*pi) - pi;
end

function plotCircle(center, radius)
  

    % Generate points for the circle
    theta = linspace(0, 2*pi, 100); % 100 points around the circle
    x = center(1) + radius * cos(theta);
    y = center(2) + radius * sin(theta);


    plot(x, y, 'k-', 'LineWidth', 1); % Circle in blue with a line width of 2
end

