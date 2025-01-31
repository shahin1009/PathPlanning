clear;close all;

% Algorithm Parameters
max_iter = 30000000;     
step_size = 15;               
goal_radius = 6; 
max_attempt = 5;

% Options
show_runtime=1; % choose 1 to show the runtime plot (slow)
profile on;

% img1 = imread("Pictures/map_1_d.png"); 
% img1 = imread("Pictures/map_2_d.png"); 
% img1 = imread("Pictures/map_3_d.png"); 
img1 = imread("Pictures/map_4_d.png"); 
% img1 = imread("Pictures/maze3.jpg");


img = imbinarize(im2gray(img1), 0.8);
se = strel('disk', 2); 
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
% goal = [271,110];

% % map2
% start= [9,79];
% goal = [262,275];
% 
% % map3
% start= [6,6];
% goal = [5,294];
% 
% % map4
% start= [11,74];
% goal = [191,178];
% 
% % map,png
% start= [65,245];
% goal = [214,148];

% maze 3
% start=[85,596];
% goal = [1056,3];


%%%%%%%%%%%%%%%%%%% Main algorithm %%%%%%%%%%%%%
tic;
path = rrt_connect_planner(img, start, goal, max_iter, step_size,max_attempt,show_runtime);
elapsedtime = toc;


if ~isempty(path)
   fprintf('Path found successfully! in %0.2f seconds\n', elapsedtime);
else
    disp('No path found.');
end

%% Bottleneck analysis
p = profile('info');
prefix = 'RRTconnectgreedy>';
function_table = p.FunctionTable;
filtered_functions = function_table(contains({function_table.FunctionName}, prefix));
for i = 1:length(filtered_functions)
    fprintf('Function: %s\n', filtered_functions(i).FunctionName);
    fprintf('Total Time: %.4f seconds\n', filtered_functions(i).TotalTime);
    fprintf('Number of Calls: %d\n', filtered_functions(i).NumCalls);
    fprintf('-----------------------------\n');
end


%%
function path = rrt_connect_planner(img, start, goal, max_iter, step_size,max_attempt,varnargin)
    % RRT Connect Planner with improved connectivity
   
    if nargin>6
        show_runtime_plot =varnargin(1);
    end

    % Initialize trees
    start_tree.nodes = start;
    start_tree.parent = 0;
    start_tree.cost = 0;
    
    goal_tree.nodes = goal;
    goal_tree.parent = 0;
    goal_tree.cost = 0;
    

    current_tree = 1;  % Alternate between trees
    
    % Main RRT Connect loop
    for iter = 1:max_iter
       
        if current_tree == 1
            q_rand = random_sampling(img, start,goal, iter, max_iter,current_tree);
            % q_rand = improved_random_sampling(img, start, goal, iter, max_iter, start_tree,goal_tree.nodes(end,:),current_tree);

            [start_tree, q_new, status1] = extend(img, start_tree, q_rand, step_size,max_attempt,show_runtime_plot);
            
            % If extension was successful
            if ~strcmp(status1, "Trapped")

                [goal_tree, connect_status] = connect(img, goal_tree, q_new, step_size,max_attempt,show_runtime_plot);
                
                % Check connection status
                if strcmp(connect_status, "Reached")
                    % Path found, reconstruct and return
                    path = reconstruct_path(start_tree, goal_tree);
                    disp(['Path found with cost: ', num2str(calculate_path_cost(path))]);
                    title(['Path found with cost: ', num2str(calculate_path_cost(path))]);
                    return;
                end
            end
            
            % Switch trees
            current_tree = 2;
        else
            q_rand = random_sampling(img, start,goal, iter, max_iter,current_tree);
            % q_rand = improved_random_sampling(img, start, goal, iter, max_iter, goal_tree,start_tree.nodes(end,:),current_tree);

            [goal_tree, q_new, status2] = extend(img, goal_tree, q_rand, step_size,max_attempt,show_runtime_plot);
            
            if ~strcmp(status2, "Trapped")
                
                [start_tree, connect_status] = connect(img, start_tree, q_new, step_size,max_attempt,show_runtime_plot);
                
                % Check connection status
                if strcmp(connect_status, "Reached")

                    path = reconstruct_path(start_tree, goal_tree);
                    disp(['Path found with cost: ', num2str(calculate_path_cost(path))]);
                    title(['Path found with cost: ', num2str(calculate_path_cost(path))]);
                    
                    return;
                end
            end
            
            % Switch trees
            current_tree = 1;
        end
    end
    
    % No path found
    path = [];
    disp('Path not found within max iterations');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%% Algorithm Helpers %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T, q_new, status] = extend(img, T, q_rand, step_size, max_attempt,show_runtime_plot)
    [nearest_idx, nearest_node] = find_nearest_node(T.nodes, q_rand);
    
    q_new = extend_node(img, nearest_node, q_rand, step_size, max_attempt);
    
    if isequal(q_new, nearest_node)
        status = "Trapped";
        return;
    end
    
 
    if norm(q_new - q_rand) <= step_size
        status = "Reached";
    else
        status = "Advanced";
    end
    

    T.nodes(end+1, :) = q_new;
    T.parent(end+1) = nearest_idx;
    if show_runtime_plot
    
        % h=plot([q_rand(1), q_new(1)], ...
        %     [q_rand(2), q_new(2)], 'r--', 'LineWidth', 0.5);
        % 
        % hh=plot(q_rand(1), q_rand(2), 'Ob', 'MarkerFaceColor', 'g');  % Circle at the start point
        % 
        % hhh=plot(q_new(1), q_new(2), 'Ob', 'MarkerFaceColor', 'g');  % Circle at the end point
        plot([nearest_node(1), q_new(1)], ...
            [nearest_node(2), q_new(2)], 'b-', 'LineWidth', 1);
    
        plot(nearest_node(1), nearest_node(2), 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
        plot(q_new(1), q_new(2), 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
        drawnow limitrate;
        % pause(1)
        % delete(h);
        % delete(hh);
        % delete(hhh);
    end
end

% attempts to connect the two trees
function [T, status] = connect(img, T, q_target, step_size,max_attempt,show_runtime_plot)

        [T, q_new, extend_status] = extend(img, T, q_target, step_size,max_attempt,show_runtime_plot);

        if strcmp(extend_status, "Reached")
            status = "Reached";
            return;
        else
        
            status = "Trapped";
            return;
        end

end

% extends a tree to a node
function new_node = extend_node(img, nearest_node, rand_point, step_size,max_attempt)

    direction = (rand_point - nearest_node) / norm(rand_point - nearest_node);
    
    new_node = nearest_node;
    
    [img_height, img_width] = size(img);
    

    
    for step = 1:max_attempt
    
        potential_new_node = new_node + direction * step_size;
        potential_new_node = round(potential_new_node);
        

        if potential_new_node(1) < 1 || potential_new_node(1) > img_width || ...
           potential_new_node(2) < 1 || potential_new_node(2) > img_height
            return;
        end
        
    
        if norm(potential_new_node - rand_point) <= step_size
           
            rand_point(1) = max(1, min(img_width, round(rand_point(1))));
            rand_point(2) = max(1, min(img_height, round(rand_point(2))));
            
            if is_path_valid(img, new_node, rand_point, step_size)
                new_node = rand_point;
                return;
            end
        end

        if ~is_path_valid(img, new_node, potential_new_node, step_size)
            return;
        end

        new_node = potential_new_node;
    end
end

function [nearest_idx, nearest_node] = find_nearest_node(nodes, point)
    % Find nearest node using Euclidean distance
    distances = sqrt(sum((nodes - point).^2, 2));
    [~, nearest_idx] = min(distances);
    nearest_node = nodes(nearest_idx, :);
end


%%%%%%%%%%%%%%%%%%%% Reconstruct the path to goal %%%%%%%%%%%%%%%%%%%%%%%%%
function path = reconstruct_path(T1, T2)
    
    conn_idx = size(T1.nodes, 1);
    
    % Reconstruct path from start
    start_path = [];
    current_idx = conn_idx;
    while current_idx > 0
        start_path = [T1.nodes(current_idx, :); start_path];
  
        if current_idx == T1.parent(current_idx)
            break;
        end
        current_idx = T1.parent(current_idx);
    end
    
    % Reconstruct path from goal
    goal_path = [];
    goal_conn_idx = size(T2.nodes, 1);
    current_idx = goal_conn_idx;
    while current_idx > 0
        goal_path = [goal_path; T2.nodes(current_idx, :)];

        if current_idx == T2.parent(current_idx)
            break;
        end
        current_idx = T2.parent(current_idx);
    end
    
    % Combine paths
    path = [start_path; goal_path];

    plot(start_path(:,1), start_path(:,2), 'r-', 'LineWidth', 3);
    plot(goal_path(:,1), goal_path(:,2), 'm-', 'LineWidth', 3);
    
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
function rand_point = improved_random_sampling(img, start, goal, iter, max_iter,T,node,current_tree)
    

    
    end_node =node;
    [~, start_node] = find_nearest_node(T.nodes, end_node);


    goal_bias_prob = 0.1;  % Probability of sampling the goal directly
    gaussian_sampling_prob = 0.2;  % Probability of Gaussian sampling

    % Random number to decide sampling strategy
    sample_strategy = rand();

    % Goal biasing: directly sample goal with a small probability
    if sample_strategy < goal_bias_prob
        if current_tree==1
            rand_point = goal;
        else
            rand_point = start;
        end
        return;
    end

    if sample_strategy < (goal_bias_prob + gaussian_sampling_prob)
        % t = rand(); % Random interpolation factor between 0 and 1
        % 
        % rand_point = round(start_node + t * (end_node - start_node));
        % 
        % % Ensure the point lies within image bounds
        % rand_point = [
        %     max(1, min(size(img, 2), rand_point(1))), ...
        %     max(1, min(size(img, 1), rand_point(2)))
        %     ];
        rand_point = end_node;

        return;
    end

    % Uniform random sampling for the rest
    rand_point = [
        randi([1, size(img,2)]), ...
        randi([1, size(img,1)])
    ];
end

function rand_point = random_sampling(img, start,goal, iter, max_iter,current_tree)

    if mod(iter,100) ==0  % Bias towards goal
        if current_tree==1
            rand_point = goal;
        else
            rand_point = start;
        end
    else
        rand_point = [randi([1, size(img,2)]), randi([1, size(img,1)])];
    end

end
