clear;
profile on;

% Algorithm Parameters
max_iter = 30000000;     
step_size = 20;        
radius = 25;          
goal_radius = 15; 

% Options
flag_show_runtime = 1;     % Showing run time plot
flag_show_randompoint = 0; % Showing random point selection

% img1 = imread("Pictures/map_1_d.png"); 
% img1 = imread("Pictures/map_2_d.png"); 
% img1 = imread("Pictures/map_3_d.png"); 
img1 = imread("Pictures/map_4_d.png"); 
% img1 = imread("Pictures/maze3.jpg");

img = imbinarize(im2gray(img1), 0.8);
se = strel('disk', 1); 
img = ~imdilate(~img, se);
map_size = size(img);

figure; % Create a figure handle
set(gcf, 'Position', [100, 100, 800, 600]); 
imshow(img, 'InitialMagnification', 'fit'); 
hold on
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
plotCircle(goal, goal_radius)

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
% 
% start=[85,596];
% goal = [1056,3];


%%%%%%%%%%%%%%%%%%% Main algorithm %%%%%%%%%%%%%
tic;
path = rrt_maze_path(img, start, goal, max_iter, step_size,radius, goal_radius,[flag_show_randompoint,flag_show_runtime]);
elapsedtime = toc;

%% Bottleneck analysis
if ~isempty(path)
   fprintf('Path found successfully! in %0.2f seconds\n', elapsedtime);
else
    disp('No path found.');
end

p = profile('info');
prefix = 'RRTOriginal>';
function_table = p.FunctionTable;
filtered_functions = function_table(contains({function_table.FunctionName}, prefix));
for i = 1:length(filtered_functions)
    fprintf('Function: %s\n', filtered_functions(i).FunctionName);
    fprintf('Total Time: %.4f seconds\n', filtered_functions(i).TotalTime);
    fprintf('Number of Calls: %d\n', filtered_functions(i).NumCalls);
    fprintf('-----------------------------\n');
end

%% Fucntion

function path = rrt_maze_path(img, start, goal, max_iter, step_size,radius, goal_radius,varnargin)
    
    
    if nargin>7
        flag_show_randompoint=varnargin(1);
        show_runtime_plot =varnargin(2);
    end
    tree.nodes = start;
    tree.parent = 0;
    tree.cost = 0;

    % Main RRT loop
    for iter = 1:max_iter
        
        % 1- Sample random point
        rand_point = random_sampling(img, goal, iter, max_iter);
        % rand_point = improved_random_sampling(img, start, goal, iter, max_iter);
        
        % 2- Find nearest node in the tree
        [nearest_idx, nearest_node] = find_nearest_node(tree.nodes, rand_point);

        % 3- Extend tree towards random point
        new_node = extend_node(nearest_node, rand_point, step_size);
        
        % 4- Check if the extention is valid
        if is_path_valid(img, nearest_node, new_node,step_size)
            
            [near_indices, near_nodes] = find_nearest_nodes(tree.nodes, new_node, radius);
            
            % Choose parent node with minimum cost path
            min_cost_idx = nearest_idx;
            min_cost = tree.cost(nearest_idx) + node_distance(nearest_node, new_node);
            
            % Check alternative near nodes for potentially lower cost path
            for i = 1:length(near_indices)
                near_idx = near_indices(i);
                near_node = near_nodes(i, :);
                
                % Check if connecting through this node is collision-free and has lower cost
                if is_path_valid(img, near_node, new_node,step_size)
                    potential_cost = tree.cost(near_idx) + node_distance(near_node, new_node);
                    if potential_cost < min_cost
                        min_cost_idx = near_idx;
                        min_cost = potential_cost;
                    end
                end
            end
            

            %5- Add new node to tree
            tree.nodes(end+1, :) = new_node;
            tree.parent(end+1) = min_cost_idx;
            tree.cost(end+1) = min_cost;


            %%%%%%%%%%%%%%%%%% visualization %%%%%%%%%%%%%%
            
            if flag_show_randompoint
                h=plot([rand_point(1), new_node(1)], ...
                    [rand_point(2), new_node(2)], 'K-', 'LineWidth', 0.5);
    
                hh=plot(rand_point(1), rand_point(2), 'ob', 'MarkerFaceColor', 'b');  % Circle at the start point
                
                hhh=plot(new_node(1), new_node(2), 'ob', 'MarkerFaceColor', 'b');  % Circle at the end point
            end
            if show_runtime_plot
                % Plot the new connection
                plot([tree.nodes(min_cost_idx,1), new_node(1)], ...
                     [tree.nodes(min_cost_idx,2), new_node(2)], 'b-', 'LineWidth', 1);
                drawnow limitrate;
            end
            if flag_show_randompoint
                delete(h);
                delete(hh);
                delete(hhh);
            end
            % Check if goal is reached
            if norm(new_node - goal) <= goal_radius
                tree.nodes(end, :) = goal;
                tree.cost(end) = tree.cost(end-1) + node_distance(goal, tree.nodes(end-1, :));
                [path,finalcost] = reconstruct_path(tree, length(tree.nodes));
                disp(['Path found with cost: ', num2str(finalcost)]);
                title(['Path found with cost: ', num2str(finalcost)]);
                return;
            end
        end
    end
    
    % No path found
    path = [];
    disp('Path not found within max iterations');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Algorithm Helpers %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nearest_idx, nearest_node] = find_nearest_node(nodes, point)
    distances = sqrt(sum((nodes - point).^2, 2));
    [~, nearest_idx] = min(distances);
    nearest_node = nodes(nearest_idx, :);
end

function [near_indices, near_nodes] = find_nearest_nodes(nodes, new_node, radius)
    distances = sqrt(sum((nodes - new_node).^2, 2));
    near_indices = find(distances <= radius);
    near_nodes = nodes(near_indices, :);
end

function new_node = extend_node(nearest_node, rand_point, step_size)
    % Extend from nearest node towards random point
    direction = (rand_point - nearest_node) / norm(rand_point - nearest_node);
    new_node = nearest_node + direction * step_size;
    new_node = round(new_node);
end


% Calculate Euclidean distance between two nodes
function dist = node_distance(node1, node2)
   
    dist = sqrt((node1(1)-node2(1))^2+((node1(2)-node2(2))^2));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% colision check %%%%%%%%%%%%%%
function is_free = is_path_valid1(img,x1, x2 ,step_size)

    num_points = round(step_size/0.4);

    
    x_points = round(linspace(round(x1(1)), round(x2(1)), num_points));
    y_points = round(linspace(round(x1(2)), round(x2(2)), num_points));

   
    x_points = max(1, min(size(img, 2), x_points)); % Clamp to [1, width]
    y_points = max(1, min(size(img, 1), y_points)); % Clamp to [1, height]

    is_free = all(arrayfun(@(x, y) img(y, x), x_points, y_points)); % 1=free, 0=obstacle
end

% Faster 
function is_free = is_path_valid(img, x1, x2, step_size)

    num_points = round(step_size/0.05);
    t = (0:num_points-1)' / (num_points-1);
    x_points = round(x1(1) + (x2(1) - x1(1)) * t);
    y_points = round(x1(2) + (x2(2) - x1(2)) * t);
    

    x_points = max(1, min(size(img, 2), x_points));
    y_points = max(1, min(size(img, 1), y_points));
    

    linear_indices = sub2ind(size(img), y_points, x_points);
    
    is_free = all(img(linear_indices));
end


function plotCircle(center, radius)

    % Generate points for the circle
    theta = linspace(0, 2*pi, 100); % 100 points around the circle
    x = center(1) + radius * cos(theta);
    y = center(2) + radius * sin(theta);


    plot(x, y, 'k-', 'LineWidth', 1); % Circle in blue with a line width of 2
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% random sampling strategies %%%%%%%%%%%%%%%%%%
function rand_point = improved_random_sampling(img, start, goal, iter, max_iter) 
    goal_bias_prob = 0.3;
    gaussian_sampling_prob = 0.2;
    
    sample_strategy = rand();
    
    if sample_strategy < goal_bias_prob
        rand_point = goal;
        return;
    end
    
    % Gaussian-like sampling around promising regions
    if sample_strategy < (goal_bias_prob + gaussian_sampling_prob)
        progress = iter / max_iter;
        range_x = max(size(img, 2) * (1 - progress), size(img, 2)/2);
        range_y = max(size(img, 1) * (1 - progress), size(img, 1)/2);
        center_x = (start(1) + goal(1)) / 2;
        center_y = (start(2) + goal(2)) / 2;
        rand_point = [
            max(1, min(size(img, 2), round(center_x + (rand()-0.5) * range_x))), ...
            max(1, min(size(img, 1), round(center_y + (rand()-0.5) * range_y)))
        ];
        
        return;
    end
    rand_point = [
        randi([1, size(img,2)]), ...
        randi([1, size(img,1)])
    ];
end

function rand_point = random_sampling(img, goal, iter, max_iter)

if mod(iter,10) ==0 
    rand_point = goal;
else
    rand_point = [randi([1, size(img,2)]), randi([1, size(img,1)])];
end
end

%%%%%%%%%%%%%%%%%%%% Reconstruct the path to goal %%%%%%%%%%%%%%%%%%%%%%%%%
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
    plot(path(:,1), path(:,2), 'r-', 'LineWidth', 3);
    
end