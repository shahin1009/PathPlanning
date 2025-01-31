clear;close all;

% Algorithm Parameters
max_iter = 10000;
step_size = 20;
radius = 25;
goal_radius = 15;
goal_bias=25;
% Options
show_runtime=1;
profile off;

img1 = imread("Pictures/map_1_d.png"); 
% img1 = imread("Pictures/map_2_d.png"); 
% img1 = imread("Pictures/map_3_d.png"); 
% img1 = imread("Pictures/map_4_d.png"); 
% img1 = imread("Pictures/maze3.jpg");

img = imbinarize(im2gray(img1), 0.8);
se = strel('disk', 1); 
img = ~imdilate(~img, se);
map_size = size(img);


figure;
set(gcf, 'Position', [100, 100, 800, 600]);
imshow(img, 'InitialMagnification', 'fit');
hold on;
title('Select START point, then GOAL point');


[start_x, start_y] = ginput(1); 
start = [round(start_x), round(start_y)];
plot(start(1), start(2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
pause(0.1); 


[goal_x, goal_y] = ginput(1);
goal = [round(goal_x), round(goal_y)];
plot(goal(1), goal(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
pause(0.1);

% map1
% start= [8,116];
% goal = [274,115];
%
% % map2
% start= [9,79];
% goal = [262,281];
% 
% % map3
% start= [6,9];
% goal = [5,298];
% 
% % map4
% start= [11,74];
% goal = [191,181];
% % 
% % map,png
% start= [65,245];
% goal = [214,148];


plotCircle(goal, goal_radius)

%%%%%%%%%%%%%%%%%%% Main algorithm %%%%%%%%%%%%%
tic;
path = rrt_star_maze_path(img, start, goal, max_iter, step_size, radius, goal_radius,[show_runtime,goal_bias]);
elapsedtime = toc;


if ~isempty(path)
    fprintf('Path found successfully! in %0.2f seconds\n', elapsedtime);
else
    disp('No path found.');
end


% profile viewer;


function path = rrt_star_maze_path(img, start, goal, max_iter, step_size, radius, goal_radius, varnargin)
if nargin > 7
    show_runtime_plot = varnargin(1);
    goal_bias = varnargin(2);
end

% Initialize tree
tree.nodes = start;
tree.parent = 0;
tree.cost = 0;

best_goal_path = [];
best_goal_cost = inf;
initial_path_found = false;
beacons = [];
n = 0;  % Iteration where initial path was found
b = goal_bias;  % Using goal_bias as the biasing ratio

% Track all nodes that reach the goal
goal_nodes_indices = [];

free_space = sum(sum(img)) / numel(img);  % Ratio of free space
total_volume = size(img,1) * size(img,2);
gamma = 2 * sqrt(free_space * total_volume / pi);

% Main RRT* loop
for iter = 1:max_iter
    if mod(iter,100)==0
        disp("iteration number:")
        disp(iter)
    end

    % Sample normally or sample around the beacons
    if initial_path_found && ~isempty(beacons) 
        if mod(iter,3)==0 
            rand_point = sample_near_beacons(beacons(2:end-1,:), radius, img);
        else
            rand_point = random_sampling(img, goal, iter, goal_bias);
        end
    else
        rand_point = random_sampling(img, goal, iter, goal_bias);
    end

    [nearest_idx, nearest_node] = find_nearest_node(tree.nodes, rand_point);
    new_node = extend_node(nearest_node, rand_point, step_size);
    current_radius = get_optimal_radius(length(tree.nodes), gamma);

    if is_path_valid(img, nearest_node, new_node, step_size)
        % Find near nodes within radius for potential connection
        [near_indices, near_nodes] = find_near_nodes(tree.nodes, new_node, current_radius);

        % Choose parent node with minimum cost path
        min_cost_idx = nearest_idx;
        min_cost = tree.cost(nearest_idx) + node_distance(nearest_node, new_node);

        % Check alternative near nodes for potentially lower cost path
        for i = 1:length(near_indices)
            near_idx = near_indices(i);
            near_node = near_nodes(i, :);

            if is_path_valid(img, near_node, new_node, step_size)
                potential_cost = tree.cost(near_idx) + node_distance(near_node, new_node);
                if potential_cost < min_cost
                    min_cost_idx = near_idx;
                    min_cost = potential_cost;
                end
            end
        end

        % Add new node to tree
        tree.nodes(end+1, :) = new_node;
        tree.parent(end+1) = min_cost_idx;
        tree.cost(end+1) = min_cost;

        if show_runtime_plot
            plot([tree.nodes(min_cost_idx,1), new_node(1)], ...
                [tree.nodes(min_cost_idx,2), new_node(2)], 'Color', [0.678, 0.847, 0.902], 'LineWidth', 1);
        end

        % Rewiring step
        for i = 1:length(near_indices)
            near_idx = near_indices(i);
            near_node = near_nodes(i, :);
            old_parent_idx = tree.parent(near_idx);

            potential_new_cost = tree.cost(end) + node_distance(new_node, near_node);

            if potential_new_cost < tree.cost(near_idx) && is_path_valid(img, new_node, near_node, step_size)
                tree.parent(near_idx) = length(tree.nodes);
                tree.cost(near_idx) = potential_new_cost;

                if show_runtime_plot
                    % Erase old connection
                    plot([tree.nodes(old_parent_idx,1), near_node(1)], ...
                        [tree.nodes(old_parent_idx,2), near_node(2)], 'color', 'w', 'LineWidth', 1.5);
                    % Draw new connection
                    plot([new_node(1), near_node(1)], ...
                        [new_node(2), near_node(2)], 'Color', [0.565, 0.933, 0.565], 'LineWidth', 1.5);
                end
            end
        end

        % Check if new node reaches goal
        if norm(new_node - goal) <= goal_radius
            goal_nodes_indices = [goal_nodes_indices, length(tree.nodes)];
            
            % Try to reconstruct path from this goal node
            [current_path, current_cost] = reconstruct_path(tree, length(tree.nodes));
            
            if ~isempty(current_path)
                if ~initial_path_found || current_cost < best_goal_cost
                    best_goal_path = current_path;
                    best_goal_cost = current_cost;
                    [beacons, path_cost_beac] = path_optimization(best_goal_path, img, step_size);
                    
                    if show_runtime_plot
                        if initial_path_found
                            % Clear previous path
                            if exist('d', 'var')
                                delete(d);
                            end
                        end
                        % Plot new best path
                        d = plot(beacons(:,1), beacons(:,2), 'm-', 'LineWidth', 3);
                        drawnow limitrate;
                    end
                    
                    if ~initial_path_found
                        initial_path_found = true;
                        n = iter;
                    end
                end
            end
        end
    end
end

% Plot and return best path
if ~isempty(best_goal_path)
    path = best_goal_path;

    if exist('d', 'var')
        delete(d);
    end
    % Plot final path
    plot(path(:,1), path(:,2), 'r-', 'LineWidth', 3);
    % Plot final beacons
    plot(beacons(:,1), beacons(:,2), 'm-', 'LineWidth', 3);

    disp(['Path found with cost: ', num2str(calculate_path_cost(beacons))]);
    title(['Path found with cost: ', num2str(calculate_path_cost(beacons))]);
else
    path = [];
    disp('Path not found within max iterations');
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%% Algorithm Helpers %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [path, actual_cost] = reconstruct_path(tree, goal_idx)
    try
        path = tree.nodes(goal_idx, :);
        current_idx = goal_idx;
        actual_cost = 0;
        visited = zeros(1, length(tree.parent));  % Track visited nodes to prevent cycles
        
        while tree.parent(current_idx) ~= 0
            % Check for cycles
            if visited(current_idx)
                path = [];
                actual_cost = inf;
                return
            end
            visited(current_idx) = 1;
            
            parent_idx = tree.parent(current_idx);
            % Add the actual segment cost
            actual_cost = actual_cost + node_distance(tree.nodes(current_idx, :), tree.nodes(parent_idx, :));
            % Add parent node to path
            current_idx = parent_idx;
            path = [tree.nodes(current_idx, :); path];
        end
    catch
        % If any error occurs during reconstruction, return empty path
        path = [];
        actual_cost = inf;
    end
end

function [optimized_path, path_cost] = path_optimization(path, img, step_size)
    % Initialize with start and goal points
    optimized_path = path(1, :);  % Start with first point
    current_idx = 1;
    path_cost = 0;
    
    while current_idx < size(path, 1)
        % Try to connect current point to furthest possible point
        for check_idx = size(path, 1):-1:current_idx + 1
            % Check if direct path exists between current point and candidate point
            if is_path_valid(img, optimized_path(end, :), path(check_idx, :), step_size)
                % Add the valid point to optimized path
                optimized_path = [optimized_path; path(check_idx, :)];
                current_idx = check_idx;
                break;
            end
            
            % If we've checked all points up to the next one
            % and found no valid path, add the next point
            if check_idx == current_idx + 1
                optimized_path = [optimized_path; path(current_idx + 1, :)];
                current_idx = current_idx + 1;
            end
        end
    end
    
    % Calculate final path cost
    for i = 1:size(optimized_path, 1)-1
        path_cost = path_cost + norm(optimized_path(i+1, :) - optimized_path(i, :));
    end
end

function rand_point = sample_near_beacons(beacons, radius, img)
% Randomly select a beacon
beacon_idx = randi(size(beacons, 1));
selected_beacon = beacons(beacon_idx, :);

% Sample within a ball of radius R_beacons centered at the beacon
angle = 2 * pi * rand();
r = radius * sqrt(rand());  % Square root for uniform distribution in circle

% Calculate point
x = round(selected_beacon(1) + r * cos(angle));
y = round(selected_beacon(2) + r * sin(angle));

% Ensure point is within image bounds
x = max(1, min(size(img, 2), x));
y = max(1, min(size(img, 1), y));

rand_point = [x, y];
end

function [nearest_idx, nearest_node] = find_nearest_node(nodes, point)
distances = sqrt(sum((nodes - point).^2, 2));
[~, nearest_idx] = min(distances);
nearest_node = nodes(nearest_idx, :);
end

function [near_indices, near_nodes] = find_near_nodes(nodes, new_node, radius)
distances = sqrt(sum((nodes - new_node).^2, 2));
near_indices = find(distances <= radius);
near_nodes = nodes(near_indices, :);
end

%steer
function new_node = extend_node(nearest_node, rand_point, step_size)
direction = (rand_point - nearest_node) / norm(rand_point - nearest_node);
new_node = nearest_node + direction * step_size;
new_node = round(new_node);
end

function dist = node_distance(node1, node2)
   
    dist = sqrt((node1(1)-node2(1))^2+((node1(2)-node2(2))^2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% colision check %%%%%%%%%%%%%%

% Faster
function is_free = is_path_valid(img, x1, x2, step_size)

num_points = round(step_size/0.01);

t = (0:num_points-1)' / (num_points-1);
x_points = round(x1(1) + (x2(1) - x1(1)) * t);
y_points = round(x1(2) + (x2(2) - x1(2)) * t);

x_points = max(1, min(size(img, 2), x_points));
y_points = max(1, min(size(img, 1), y_points));

linear_indices = sub2ind(size(img), y_points, x_points);

is_free = all(img(linear_indices));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%% random sampling strategies %%%%%%%%%%%%%%%%%%

function rand_point = random_sampling(img, goal, iter, goal_bias)

    if nargin < 3
        goal_bias =  75;
    end
    
    
    if mod(iter,goal_bias) ==0  % Bias towards goal
        rand_point = goal;
    else
        rand_point = [randi([1, size(img,2)]), randi([1, size(img,1)])];
    end

end

%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%

function radius = get_optimal_radius(n, gamma)

if n < 2
    radius = inf; 
    return;
end

d = 2;  % 2D space
radius = gamma * (log(n)/n)^(1/d);
end


function plotCircle(center, radius)

% Validate inputs
if length(center) ~= 2 || ~isnumeric(center)
    error('Center must be a numeric vector of size 1x2.');
end
if ~isscalar(radius) || ~isnumeric(radius) || radius <= 0
    error('Radius must be a positive numeric scalar.');
end

theta = linspace(0, 2*pi, 100); % 100 points around the circle
x = center(1) + radius * cos(theta);
y = center(2) + radius * sin(theta);


plot(x, y, 'k-', 'LineWidth', 1); % Circle in blue with a line width of 2
end
