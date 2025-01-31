clear; close all;


% Main function to run PRM path planning

% Load and prepare the map
img1 = imread("Pictures/map_3_d.png");

% map1
% start= [8,116];
% goal = [274,115];

% % map2
% start= [9,79];
% goal = [262,281];

%
% % map3
start= [6,6];
goal = [5,298];
%
% % map4
% start= [11,74];
% goal = [191,181];
%
% % map,png
% start= [65,245];
% goal = [214,148];


img = imbinarize(im2gray(img1), 0.8);
se = strel('disk', 1);
img = ~imdilate(~img, se);
map_size = size(img);

% Get parameters and map
params = initialize_parameters();

figure;
set(gcf, 'Position', [100, 100, 800, 600]);
imshow(img, 'InitialMagnification', 'fit');
hold on
title('PRM Path Planning');

% Generate PRM
nodes = generate_prm_nodes(img, map_size, params.n_nodes, start, goal);

% Build roadmap
adjacency_matrix = build_roadmap(nodes, img, map_size, params.k_nearest);

start_node = size(nodes, 1) - 1;
goal_node = size(nodes, 1);
% Find path (Dijkstra)
path = find_path(nodes, adjacency_matrix, start_node, goal_node);

for i = 1:(length(path)-1)
        plot([nodes(path(i),1), nodes(path(i+1),1)], ...
            [nodes(path(i),2), nodes(path(i+1),2)], ...
            'Color', 'r', 'LineWidth', 2);
end



function params = initialize_parameters()
    % Initialize algorithm parameters
    params.flag_show_runtime = 1;
    params.n_nodes = 1000;
    params.k_nearest = 20;
end


function nodes = generate_prm_nodes(img, map_size, n_nodes, start, goal)
    % Generate random nodes for PRM
    nodes = zeros(n_nodes, 2);
    valid_nodes = 0;
    
    while valid_nodes < n_nodes
        x = randi(map_size(2));
        y = randi(map_size(1));
        
        if img(y, x) == 1
            valid_nodes = valid_nodes + 1;
            nodes(valid_nodes, :) = [x, y];
            plot(x, y, 'g.', 'MarkerSize', 5);
        end
    end
    
    % Add start and goal to nodes
    nodes = [nodes; start(1:2); goal(1:2)];
    plot(start(1), start(2), 'ro', 'MarkerSize', 10,'MarkerFaceColor', 'r');
    plot(goal(1), goal(2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
end

function adjacency_matrix = build_roadmap(nodes, img, map_size, k_nearest)
    % Build the roadmap by connecting nodes
    n_total = size(nodes, 1);
    adjacency_matrix = zeros(n_total);
    
    for i = 1:n_total
        distances = sqrt(sum((nodes - nodes(i,:)).^2, 2));
        [~, nearest_idx] = sort(distances);
        
        for j = 2:(k_nearest+1)
            if j > length(nearest_idx)
                break;
            end
            neighbor_idx = nearest_idx(j);
            
            if is_path_free(nodes(i,:), nodes(neighbor_idx,:), img, map_size)
                adjacency_matrix(i, neighbor_idx) = 1;
                adjacency_matrix(neighbor_idx, i) = 1;
                
                node1=nodes(i,:);
                node2=nodes(neighbor_idx,:);
                %Plot
                line([node1(1), node2(1)], [node1(2), node2(2)], 'Color', [0.7 0.7 0.7]);
                % drawnow limitrate;
            end
        end
    end
end

function free = is_path_free(node1, node2, img, map_size)


    path_x = linspace(node1(1), node2(1), 20);
    path_y = linspace(node1(2), node2(2), 20);
    
    path_x_checked = min(max(round(path_x), 1), map_size(2));
    path_y_checked = min(max(round(path_y), 1), map_size(1));
    
    free = true;
    for k = 1:length(path_x)
        if img(path_y_checked(k), path_x_checked(k)) == 0
            free = false;
            break;
        end
    end
end


% Dijkstra
function path = find_path(nodes, adjacency_matrix, start_idx, goal_idx)
    % Find path using Dijkstra's algorithm
    n_total = size(nodes, 1);
    distances = inf(n_total, 1);
    distances(start_idx) = 0;
    parent = zeros(n_total, 1);
    visited = false(n_total, 1);
    
    while ~all(visited)
        % New candidate
        unvisited_distances = distances;
        unvisited_distances(visited) = inf;
        [~, current] = min(unvisited_distances);
        
        if current == goal_idx
            break;
        end
        
        neighbors = find(adjacency_matrix(current, :));
        for neighbor = neighbors
            if ~visited(neighbor)
                distance = distances(current) + ...
                    norm(nodes(current,:) - nodes(neighbor,:));
                if distance < distances(neighbor)
                    % Update Info
                    distances(neighbor) = distance;
                    parent(neighbor) = current;
                end
            end
        end
        
        visited(current) = true;
    end
    
    % Reconstruct path
    path = reconstruct_path(parent, start_idx, goal_idx);
end

function path = reconstruct_path(parents, start_idx, goal_idx)
    % Reconstruct path from parents array
    path = goal_idx;
    current = goal_idx;
    while current ~= start_idx
        current = parents(current);
        path = [current, path];
    end
end



