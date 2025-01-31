clear;close all;

%%
% Algorithm Parameters
max_iter = 20000;     
step_size = 15;        
radius = 10;          
goal_radius = 10; 
goal_bias=20;
% Options
show_runtime=1;
profile off;

% img1 = imread("Pictures/map_1_d.png"); 
img1 = imread("Pictures/map_2_d.png"); 
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


% select by the user
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

% % map2
% start= [9,79];
% goal = [262,281];
% 
% % map3
% start= [6,6];
% goal = [5,298];
% 
% % map4
% start= [11,74];
% goal = [191,181];
% 
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
    if nargin>7
        show_runtime_plot = varnargin(1);
        goal_bias = varnargin(2);
    end

  
    d=[];
    % Initialize tree
    tree.nodes = start;
    tree.parent = 0;
    tree.cost = 0;
    initial_path_found = false;
    best_goal_path = [];
    best_cost = inf;
    c_best = inf;  % Best cost so far
    
    % Calculate the distance between start and goal
    c_min = node_distance(start, goal);

    free_space = sum(sum(img)) / numel(img);  % Ratio of free space
    total_volume = size(img,1) * size(img,2);
    gamma = 2 * sqrt(free_space * total_volume / pi); 
    
    % Main Informed RRT* loop
    for iter = 1:max_iter

        if mod(iter,100)==0
            disp("iteration number:")
            disp(iter)
        end
        % If we have a solution, sample from the informed subset
        if c_best < inf
            rand_point = informed_sample(start, goal, c_best, c_min, img);
        else
            rand_point = random_sampling(img, goal, iter, goal_bias);
        end

        [nearest_idx, nearest_node] = find_nearest_node(tree.nodes, rand_point);
        new_node = extend_node(nearest_node, rand_point, step_size);
        current_radius = get_optimal_radius(length(tree.nodes), gamma);
        current_radius = max(current_radius,20);
        if is_path_valid(img, nearest_node, new_node, step_size)
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
                        plot([tree.nodes(old_parent_idx,1), near_node(1)], ...
                            [tree.nodes(old_parent_idx,2), near_node(2)], 'color', 'w', 'LineWidth', 1.5);
                        plot([new_node(1), near_node(1)], ...
                            [new_node(2), near_node(2)], 'Color', [0.565, 0.933, 0.565], 'LineWidth', 1.5);


                    end
                end
            end
            
            % Check if goal is reached
            if norm(new_node - goal) <= goal_radius && ~initial_path_found
                goal_ind = length(tree.nodes);
                current_path_cost = tree.cost(end) + node_distance(new_node, goal);
                initial_path_found = true;
                % Update best path if this path is better
                if current_path_cost < c_best
                    c_best = current_path_cost;
                    tree.nodes(end, :) = goal;
                    tree.cost(end) = current_path_cost;
                    [best_goal_path,best_goal_cost] = reconstruct_path(tree, length(tree.nodes));
                    if ~isempty(d)
                        delete(d);
                    end
                    d=plot(best_goal_path(:,1), best_goal_path(:,2), 'r-', 'LineWidth', 3);
                    
                    if show_runtime_plot
                        % Plot the ellipse
                        plot_informed_set(start, goal, c_best, c_min);
                      
                    end
                end
            elseif initial_path_found
                [current_path,current_cost] = reconstruct_path(tree, goal_ind);
                best_goal_path=current_path;
                if current_cost<best_goal_cost
                    best_goal_cost=current_cost;
                    c_best = current_cost;
                    if ~isempty(d)
                        delete(d);
                    end
                    d=plot(best_goal_path(:,1), best_goal_path(:,2), 'r-', 'LineWidth', 3);
                    
                    if show_runtime_plot
                        % Plot the ellipse
                        plot_informed_set(start, goal, c_best, c_min);
                     
                    end
                end
            end
        end
    end
    
    % Plot and return best path
    if ~isempty(best_goal_path)
        delete(d);
        path = best_goal_path;
        plot(path(:,1), path(:,2), 'r-', 'LineWidth', 3);
        drawnow;
        disp(['Path found with cost: ', num2str(c_best)]);
        title(['Path found with cost: ', num2str(best_goal_cost)]);
    else
        path = [];
        disp('Path not found within max iterations');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%% Algorithm Helpers %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_informed_set(start, goal, c_best, c_min)

    center = (start + goal) / 2;
    direction = (goal - start) / norm(goal - start);
    angle = atan2(direction(2), direction(1));
    

    a = c_best / 2;
    b = sqrt(c_best^2 - c_min^2) / 2;
    

    theta = linspace(0, 2*pi, 100);
    x = a * cos(theta);
    y = b * sin(theta);

    rotation = [cos(angle) -sin(angle); sin(angle) cos(angle)];
    points = rotation * [x; y];
    

    x_ellipse = points(1, :) + center(1);
    y_ellipse = points(2, :) + center(2);
    
    plot(x_ellipse, y_ellipse, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    drawnow;
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
    % Calculate number of points
    num_points = round(step_size/0.1);
    
 
    t = (0:num_points-1)' / (num_points-1);
    x_points = round(x1(1) + (x2(1) - x1(1)) * t);
    y_points = round(x1(2) + (x2(2) - x1(2)) * t);
    

    x_points = max(1, min(size(img, 2), x_points));
    y_points = max(1, min(size(img, 1), y_points));
    
    linear_indices = sub2ind(size(img), y_points, x_points);
  
    is_free = all(img(linear_indices));
end

%%%%%%%%%%%%%%%%%%%% Reconstruct the path to goal %%%%%%%%%%%%%%%%%%%%%%%%%

function path = reconstruct_path2(tree, goal_idx)
    
    path = tree.nodes(goal_idx, :);
    current_idx = goal_idx;
    
    while tree.parent(current_idx) ~= 0
        current_idx = tree.parent(current_idx);
        path = [tree.nodes(current_idx, :); path];
    end
    
    
end

function [path, actual_cost] = reconstruct_path(tree, goal_idx)
    path = tree.nodes(goal_idx, :);
    current_idx = goal_idx;
    actual_cost = 0;
    
    while tree.parent(current_idx) ~= 0
        parent_idx = tree.parent(current_idx);
        % Add the actual segment cost
        actual_cost = actual_cost + node_distance(tree.nodes(current_idx, :), tree.nodes(parent_idx, :));
        % Add parent node to path
        current_idx = parent_idx;
        path = [tree.nodes(current_idx, :); path];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% random sampling strategies %%%%%%%%%%%%%%%%%%

function rand_point = random_sampling(img, goal, iter, goal_bias)

if nargin < 3
    goal_bias =  10;
end


if mod(iter,goal_bias) ==0  % Bias towards goal
    rand_point = goal;
else
    rand_point = [randi([1, size(img,2)]), randi([1, size(img,1)])];
end

end

function sample = informed_sample(start, goal, c_best, c_min, img)
    if abs(c_best - c_min) < 0.1
        sample = random_sampling(img, goal, 1, 20);
        return;
    end
    
    % Calculate center point
    center = (start + goal) / 2;
    
    % Calculate rotation matrix to align ellipse
    direction = (goal - start) / norm(goal - start);
    angle = atan2(direction(2), direction(1));
    rotation = [cos(angle) -sin(angle); sin(angle) cos(angle)];
    
    % Calculate ellipse dimensions
    a = c_best / 2;
    b = sqrt(c_best^2 - c_min^2) / 2;
    
   
    while true
        
        theta = 2 * pi * rand();
        r = rand();
        
        x = a * r * cos(theta);
        y = b * r * sin(theta);
        
        % Rotate and translate
        point = rotation * [x; y] + center';
        
        sample = round(point');
        if sample(1) >= 1 && sample(1) <= size(img, 2) && ...
           sample(2) >= 1 && sample(2) <= size(img, 1)
            break;
        end
    end
end

%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%

function radius = get_optimal_radius(n, gamma)

    if n < 2
        radius = inf;  % Use maximum radius for first few nodes
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

    % Generate points for the circle
    theta = linspace(0, 2*pi, 100); % 100 points around the circle
    x = center(1) + radius * cos(theta);
    y = center(2) + radius * sin(theta);


    plot(x, y, 'k-', 'LineWidth', 1); % Circle in blue with a line width of 2
end


function plotelips()

theta = linspace(0, 2*pi, 100);
x = 2 * cos(theta);
y = 5 * sin(theta);
angle=pi/3;
center=[0,0];

rotation = [cos(angle) -sin(angle); sin(angle) cos(angle)];
points = rotation * [x; y];

x_ellipse = points(1, :) + center(1);
y_ellipse = points(2, :) + center(2);

plot(x_ellipse, y_ellipse, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);hold on;
plot(x, y, '--', 'Color', 'r', 'LineWidth', 1);
drawnow;
end