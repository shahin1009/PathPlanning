clear;close all;

% Algorithm Parameters
max_iter = 30000000;
step_size = 20;
goal_radius = 15;
v_max = 4;     
omega_max = 8; 
min_turning_radius = 23.0;  

profile off;
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
% img1 = imread("Pictures/map.pgm");
% start= [65,245,3*pi/2];
% goal = [214,148,pi/2];

% maze 3
img1 = imread("Pictures/maze3.jpg");
start=[85,596,3*pi/2];
goal = [1056,3,3*pi/2];


img = imbinarize(im2gray(img1), 0.8);
se = strel('disk', 1); 
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
path = blossom_kinodynamic_rrt(img, start, goal, max_iter, step_size, goal_radius,v_max,omega_max,min_turning_radius,flag_show_runtime);
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


function path = rrt_blossom(img, start, goal, max_iter, step_size, goal_radius, v_max, omega_max, min_turning_radius, flag_show_runtime)
    % Initialize tree structure
    tree.nodes = start;
    tree.parent = 0;
    tree.cost = 0;
    tree.controls = [];
    tree.trajectories = cell(1);
    tree.path_costs = 0;
    tree.viability = {'untried'}; % Track viability status of edges
    
    % Plot start and goal
    plot(start(1), start(2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(goal(1), goal(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    
    % Initialize parameters
    robot_radius = 1;
    [img_height, img_width] = size(img);
    dormant_deadlock = false;
    
    for iter = 1:max_iter
        if mod(iter,1000)==0
            disp(iter);
        end
        
        % Sample random point with goal bias
        if rand() < 0.1
            rand_point = goal;
        else
            rand_point = [randi([1, img_width]), randi([1, img_height]), (rand() * 2) * pi];
        end
        
        % Find nearest node
        [nearest_idx, nearest_node] = find_nearest_node(tree.nodes, rand_point);
        
        % Node blossom - try all possible controls
        new_nodes = [];
        new_trajectories = {};
        new_controls = [];
        new_costs = [];
        
        for i = 1:8  % Try multiple control samples
            control = sample_control(v_max, omega_max, min_turning_radius);
            [new_state, trajectory] = generate_trajectory_euler(nearest_node, control, step_size);
            
            % Check if trajectory is collision-free
            if ~is_trajectory_collision_free(img, trajectory, robot_radius)
                continue;
            end
            
            % Check for regression unless in deadlock
            if ~dormant_deadlock && is_regressing(tree, nearest_idx, new_state)
                continue;
            end
            
            % Store valid node
            new_nodes = [new_nodes; new_state];
            new_trajectories{end+1} = trajectory;
            new_controls = [new_controls; control];
            new_costs = [new_costs; compute_trajectory_cost(trajectory)];
        end
        
        % Add all valid nodes to tree
        for i = 1:size(new_nodes, 1)
            current_idx = size(tree.nodes, 1) + 1;
            tree.nodes(current_idx, :) = new_nodes(i,:);
            tree.parent(current_idx) = nearest_idx;
            tree.controls(current_idx, :) = new_controls(i,:);
            tree.trajectories{current_idx} = new_trajectories{i};
            tree.path_costs(current_idx) = tree.path_costs(nearest_idx) + new_costs(i);
            tree.viability{current_idx} = 'live';
            
            % Plot the trajectory
            if flag_show_runtime
                plot(new_trajectories{i}(:,1), new_trajectories{i}(:,2), 'b-', 'LineWidth', 0.5);
                drawnow limitrate;
            end
            
            % Check if goal reached
            if norm(new_nodes(i,1:2) - goal(1:2)) <= goal_radius
                [path, final_cost] = reconstruct_path(tree);
                fprintf('Final path cost: %.2f\n', final_cost);
                return;
            end
        end
        
        % Check for dormant deadlock
        if is_dormant_deadlock(tree)
            dormant_deadlock = true;
        end
    end
    
    path = [];
    warning('Path not found within maximum iterations');
end

function is_regressing = is_regressing(tree, parent_idx, new_state)
    % Check if new state regresses into already explored space
    is_regressing = false;
    parent_node = tree.nodes(parent_idx, :);
    
    for i = 1:size(tree.nodes, 1)
        if i == parent_idx
            continue;
        end
        % Only consider live nodes for regression check
        if strcmp(tree.viability{i}, 'live')
            if compute_distance(tree.nodes(i,:), new_state) < ...
               compute_distance(parent_node, new_state)
                is_regressing = true;
                return;
            end
        end
    end
end

function deadlock = is_dormant_deadlock(tree)
    % Check if all branches are either dead or dormant
    deadlock = true;
    for i = 1:length(tree.viability)
        if strcmp(tree.viability{i}, 'live') || strcmp(tree.viability{i}, 'untried')
            deadlock = false;
            return;
        end
    end
end

function update_viability(tree, node_idx)
    % Update viability status based on node's children
    children_idx = find(tree.parent == node_idx);
    all_children_dead = true;
    all_children_dormant = true;
    
    for child_idx = children_idx
        status = tree.viability{child_idx};
        if ~strcmp(status, 'dead')
            all_children_dead = false;
        end
        if ~strcmp(status, 'dormant')
            all_children_dormant = false;
        end
    end
    
    if all_children_dead
        tree.viability{node_idx} = 'dead';
    elseif all_children_dormant
        tree.viability{node_idx} = 'dormant';
    end
end

% The following helper functions remain unchanged:
% - generate_trajectory_euler
% - compute_trajectory_cost
% - compute_distance
% - find_nearest_node
% - sample_control
% - is_trajectory_collision_free
% - reconstruct_path
% - angle_diff

function path = blossom_kinodynamic_rrt(img, start, goal, max_iter, step_size, goal_radius, v_max, omega_max, min_turning_radius, flag_show_runtime)
    % Initialize tree structure
    tree.nodes = start;    % [x, y, theta]
    tree.parent = 0;       % Parent index
    tree.cost = 0;         % Cost from start
    tree.controls = [];    % Store control inputs
    tree.trajectories = cell(1); % Store complete trajectories
    tree.path_costs = 0;   % Add path_costs to store distance-based costs

    % Blossom parameters
    num_blossoms = 3;      % Number of candidate paths to generate at each iteration
    num_control_samples = 5; % Number of control samples for each blossom
    
    % Plot start and goal
    plot(start(1), start(2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(goal(1), goal(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    
    % Initialize collision checking parameters
    robot_radius = 2;
    [img_height, img_width] = size(img);
    
    for iter = 1:max_iter
        % 1. Sample multiple random points (blossoms)
        blossom_points = generate_blossom_points(num_blossoms, img_width, img_height, goal);
        
        if mod(iter,1000)==0
            disp(iter);
        end

        % For each blossom point, try to grow the tree
        for b = 1:num_blossoms
            rand_point = blossom_points(b, :);
            
            % 2. Find nearest node
            [nearest_idx, nearest_node] = find_nearest_node(tree.nodes, rand_point);
            
            % 3. Generate multiple controls and select best
            best_state = [];
            best_control = [];
            best_trajectory = [];
            best_path_cost = inf;
            min_dist = inf;
            
            % Try multiple control samples for each blossom
            for i = 1:num_control_samples
                control = sample_control(v_max, omega_max, min_turning_radius);
                [new_state, trajectory] = generate_trajectory_euler(nearest_node, control, step_size);
                
                if is_trajectory_collision_free(img, trajectory, robot_radius)
                    % Compute weighted cost combining distance and smoothness
                    dist = compute_distance(new_state, rand_point);
                    path_cost = compute_trajectory_cost(trajectory);
                    smoothness_cost = compute_smoothness_cost(trajectory);
                    total_cost = dist + 0.3 * path_cost + 0.2 * smoothness_cost;
                    
                    if total_cost < min_dist
                        min_dist = total_cost;
                        best_state = new_state;
                        best_control = control;
                        best_trajectory = trajectory;
                        best_path_cost = path_cost;
                    end
                end
            end
            
            % Add the best candidate to the tree
            if ~isempty(best_state)
                current_idx = size(tree.nodes, 1) + 1;
                tree.nodes(current_idx, :) = best_state;
                tree.parent(current_idx) = nearest_idx;
                tree.controls(current_idx, :) = best_control;
                tree.path_costs(current_idx) = tree.path_costs(nearest_idx) + best_path_cost;
                tree.trajectories{current_idx} = best_trajectory;
                
                % Plot the trajectory
                if flag_show_runtime
                    plot(best_trajectory(:,1), best_trajectory(:,2), 'b-', 'LineWidth', 0.5);
                    drawnow limitrate;
                end
                
                % Check if goal reached
                if norm(best_state(1:2) - goal(1:2)) <= goal_radius
                    [path, final_cost] = reconstruct_path(tree);
                    fprintf('Final path cost: %.2f\n', final_cost);
                    return;
                end
            end
        end
    end
    
    path = [];
    warning('Path not found within maximum iterations');
end

function blossom_points = generate_blossom_points(num_points, img_width, img_height, goal)
    blossom_points = zeros(num_points, 3);
    
    % Goal biasing probability
    goal_bias = 0.2;
    
    for i = 1:num_points
        if rand() < goal_bias
            blossom_points(i, :) = goal;
        else
            % Generate random point with Gaussian distribution around the current region
            x = randi([1, img_width]);
            y = randi([1, img_height]);
            theta = (rand() * 2) * pi;
            blossom_points(i, :) = [x, y, theta];
        end
    end
end

function smoothness_cost = compute_smoothness_cost(trajectory)
    smoothness_cost = 0;
    
    % Compute smoothness based on change in heading
    for i = 2:size(trajectory,1)-1
        % Calculate change in heading between consecutive segments
        v1 = trajectory(i,:) - trajectory(i-1,:);
        v2 = trajectory(i+1,:) - trajectory(i,:);
        
        % Calculate angle between vectors (ignoring theta component)
        angle = atan2(v2(2), v2(1)) - atan2(v1(2), v1(1));
        angle = abs(angle_diff(angle, 0));
        
        smoothness_cost = smoothness_cost + angle;
    end
end


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

    if abs(omega)<0.06
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

