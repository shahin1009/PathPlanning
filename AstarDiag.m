% Astar Diagonal with bigger with sparse edge matrix
clear;
close all;


% img1 = imread("Pictures/map_1_d.png"); 
% img1 = imread("Pictures/map_2_d.png"); 
img1 = imread("Pictures/map_3_d.png"); 
% img1 = imread("Pictures/map_4_d.png"); 
% img1 = imread("Pictures/maze3.jpg");

scaleFactor = 1; % Increase resolution by a factor of 2
img = imresize(img1, scaleFactor, 'bilinear'); % or 'nearest' for binary images

img = imbinarize(im2gray(img), 0.8);
se = strel('disk', 2); 
img = ~imdilate(~img, se);
map_size = size(img);
plotruntime = 0; % choose 1 to show the runtime plot (slow)
% Heuristic function choice
func = @ heuristic1;  
alfa=0.3; % lower: more greedy

figure; % Create a figure handle
set(gcf, 'Position', [100, 100, 800, 600]); 
imshow(img, 'InitialMagnification', 'fit'); 
hold on
title('Select START point, then GOAL point');


[start_x, start_y] = ginput(1); 
start = [round(start_x), round(start_y)];
plot(start(1), start(2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
pause(0.1); 

[goal_x, goal_y] = ginput(1);
goal = [round(goal_x), round(goal_y)];
plot(goal(1), goal(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
pause(0.1);
profile on;

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

plot(start(1), start(2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(goal(1), goal(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

%%
pathw1 = Astarcode(start,goal,img,plotruntime,func,alfa,scaleFactor);

%% Bottleneck analysis
% 
% p = profile('info');
% prefix = 'RRTOriginal>';
% function_table = p.FunctionTable;
% filtered_functions = function_table(contains({function_table.FunctionName}, prefix));
% for i = 1:length(filtered_functions)
%     fprintf('Function: %s\n', filtered_functions(i).FunctionName);
%     fprintf('Total Time: %.4f seconds\n', filtered_functions(i).TotalTime);
%     fprintf('Number of Calls: %d\n', filtered_functions(i).NumCalls);
%     fprintf('-----------------------------\n');
% end

%% Functions
function pathw = Astarcode(start_pos,goal_pos,BW,plotruntime,func,alfa,scaleFactor)

% initialize variables
dist=inf*ones(size(BW,1)*size(BW,2),1);
prec=inf*ones(size(BW,1)*size(BW,2),1); 
nodelist=-1*ones(size(BW,1)*size(BW,2),1); % to visit
f_score = inf*ones(size(BW,1)*size(BW,2),1); 
G = createSparseGraph(BW,scaleFactor);



start= (start_pos(2)-1)*size(BW,2)+start_pos(1);
goal= (goal_pos(2)-1)*size(BW,2)+goal_pos(1);

%%% visualize map start & goal
plot(start_pos(1),start_pos(2),'sg','MarkerFaceColor','g')
plot(goal_pos(1),goal_pos(2),'sr','MarkerFaceColor','r')



%%%% A* Algorithm %%%%%%%%%%%%%%%%%%%%
%%% initialize actual node = start
dist(start) = 0;
f_score(start) = func(start, goal, BW,scaleFactor);
act_node = start;
[~,con_nodes] = find(G(act_node,:)>0);
nodelist(con_nodes,1) = 1;
nodelist(act_node,1) = 0;
tic;
%%% Main Loop
while any(nodelist(:,1)==1) && act_node~=goal
    i_con = length(con_nodes);
    while i_con>0 
        edge_cost = G(act_node, con_nodes(i_con));  % Get actual cost from G matrix
        tentative_g_score = dist(act_node) + edge_cost;
        if tentative_g_score < dist(con_nodes(i_con))
            dist(con_nodes(i_con)) = tentative_g_score;
            prec(con_nodes(i_con)) = act_node;
            % f_score(con_nodes(i_con)) = alfa*dist(con_nodes(i_con)) + (1-alfa)*func(con_nodes(i_con), goal, BW,scaleFactor);
            f_score(con_nodes(i_con)) = dist(con_nodes(i_con)) + func(con_nodes(i_con), goal, BW,scaleFactor);

        end
        i_con = i_con-1;
    end

    %%% plot run-time
    if plotruntime
        plot_runtime(act_node,BW)
    end

    %%% evaluate new candidate node & new neighbours
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    act_node = new_candidate(nodelist,f_score,goal,BW,func,scaleFactor);
    % act_node = new_candidate2(nodelist,f_score,goal,BW);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nodelist(act_node,1) = 0; % visited
    [~,con_nodes] = find(G(act_node,:)>0);
    i_con = length(con_nodes);
    while i_con>0
        if nodelist(con_nodes(i_con),1) ~= 0 % if not visited
            nodelist(con_nodes(i_con),1) = 1; % add to visit
        end
        i_con = i_con-1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%% Results %%%%%%%%%%%%%%%%%%%%%%%%%%

elapsedTime = toc; % Stop timer and get elapsed time
disp(['Elapsed time: ', num2str(elapsedTime), ' seconds']);

%%% shortest path (same as before)
if dist(goal)<inf
    sol_id = goal;
    path = [];
    while(sol_id~=start)
        path = [sol_id path];
        sol_id = prec(sol_id);
    end
    path = [sol_id path];
    %%% plot shortest path
    pathw=zeros(length(path),2);
    for i=1:length(path)
        pathw(i,:)=[mod(path(i),size(BW,2)), floor(path(i)/size(BW,1))+1];
        plot(mod(path(i),size(BW,2)), floor(path(i)/size(BW,1))+1,'or','MarkerFaceColor','r');
        
    end
    title(['Total path cost: ', num2str(dist(goal))]);
    fprintf('Total path cost: %.2f\n', dist(goal));
    fprintf('Total number of nodes: %.2f\n', size(path,2));

else
    pathw=[];
    disp('no solution found')
end

disp('total number of explored nodes:')
disp(numel(nodelist) - nnz(nodelist))


end
%% Functions
%%%%%%%%%%%%%%%%%% Distance Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = heuristic1(node, goal, BW,scaleFactor)
    [node_i, node_j] = ind2sub(size(BW), node);
    [goal_i, goal_j] = ind2sub(size(BW), goal);
    h = (abs(node_i - goal_i) + abs(node_j - goal_j))/ scaleFactor;
end

function h = heuristic2(node, goal, BW,scaleFactor)
    [node_i, node_j] = ind2sub(size(BW), node);
    [goal_i, goal_j] = ind2sub(size(BW), goal);
    h = (sqrt((node_i - goal_i)^2 + (node_j - goal_j)^2))/ scaleFactor;
end

function h = heuristic3(node, goal, BW,scaleFactor)
    [node_i, node_j] = ind2sub(size(BW), node);
    [goal_i, goal_j] = ind2sub(size(BW), goal);
    h = max(abs(node_i - goal_i), abs(node_j - goal_j))/ scaleFactor;
end

function h = heuristic4(node, goal, BW,scaleFactor)
    [node_i, node_j] = ind2sub(size(BW), node);
    [goal_i, goal_j] = ind2sub(size(BW), goal);
    dx = abs(node_i - goal_i);
    dy = abs(node_j - goal_j);
    h = max(dx, dy) + (sqrt(2) - 1) * min(dx, dy)/ scaleFactor;
end

%%%%%%%%%%%%% New candidate selection %%%%%%%%%%%%%%%%%%%%%%
function act_node = new_candidate(nodelist,f_score,goal,BW,func,scaleFactor)

valid_nodes = find(nodelist(:,1)==1);

valid_f_scores = f_score(valid_nodes);
min_f_score = min(valid_f_scores);
min_f_nodes = valid_nodes(valid_f_scores == min_f_score);

if length(min_f_nodes) > 1
    % Calculate actual distances to goal for tied nodes
    distances_to_goal = zeros(size(min_f_nodes));
    for i = 1:length(min_f_nodes)
        node = min_f_nodes(i);
        
        distances_to_goal(i) =  func(node, goal, BW,scaleFactor);
    end
    [~, min_idx] = min(distances_to_goal);
    act_node = min_f_nodes(min_idx);
else
    act_node = min_f_nodes(1);
end

end


function act_node = new_candidate2(nodelist,f_score,goal,BW)
    [min_val,~] = min(f_score(nodelist(:,1)==1));
    new_nodes = find(f_score==min_val);
    tmp_i = 1;
    while nodelist(new_nodes(tmp_i),1)~=1
        tmp_i = tmp_i+1; % first best
    end
    act_node = new_nodes(tmp_i); % new node
end

%%%%%%%%%%%%%%%%% Edge matrix creation %%%%%%%%%%%%%%%%%%%%%%
function G = createSparseGraph2(BW)
    % Find indices of non-zero elements in BW
    [rows, cols] = find(BW);
    
    % Preallocate sparse matrix indices - 8 possible neighbors per node
    I = zeros(length(rows)*8, 1);
    J = zeros(length(rows)*8, 1);
    V = zeros(length(rows)*8, 1);
    
    count = 0;
    for k = 1:length(rows)
        i = rows(k);
        j = cols(k);
        
        % Define all 8 possible neighbors
        neighbors = [
            i+1, j;    % below
            i-1, j;    % above
            i, j+1;    % right
            i, j-1;    % left
            i-1, j-1;  % diagonal up-left
            i-1, j+1;  % diagonal up-right
            i+1, j-1;  % diagonal down-left
            i+1, j+1   % diagonal down-right
        ];
        
        % Check each neighbor
        for n = 1:size(neighbors, 1)
            ni = neighbors(n,1);  % neighbor row
            nj = neighbors(n,2);  % neighbor column
            
            % Check if neighbor is within bounds and is traversable (=1 in BW)
            if ni > 0 && ni <= size(BW,1) && ...
               nj > 0 && nj <= size(BW,2) && ...
               BW(ni,nj) == 1
                
                count = count + 1;
                % Convert 2D indices to 1D for sparse matrix
                I(count) = (i-1)*size(BW,2) + j;        % current node
                J(count) = (ni-1)*size(BW,2) + nj;      % neighbor node
                
                % Calculate true Euclidean distance
                dx = abs(ni - i);
                dy = abs(nj - j);
                V(count) = sqrt(dx^2 + dy^2);
            end
        end
    end
    
    % Trim excess preallocated space
    I = I(1:count);
    J = J(1:count);
    V = V(1:count);
    
    % Create sparse matrix
    G = sparse(I, J, V, size(BW,1)*size(BW,2), size(BW,1)*size(BW,2));
end

function G = createSparseGraph(BW, scaleFactor)
    % Find indices of non-zero elements in BW
    [rows, cols] = find(BW);
    
    % Preallocate sparse matrix indices - 8 possible neighbors per node
    I = zeros(length(rows)*8, 1);
    J = zeros(length(rows)*8, 1);
    V = zeros(length(rows)*8, 1);
    
    count = 0;
    for k = 1:length(rows)
        i = rows(k);
        j = cols(k);
        
        % Define all 8 possible neighbors
        neighbors = [
            i+1, j;    % below
            i-1, j;    % above
            i, j+1;    % right
            i, j-1;    % left
            i-1, j-1;  % diagonal up-left
            i-1, j+1;  % diagonal up-right
            i+1, j-1;  % diagonal down-left
            i+1, j+1   % diagonal down-right
        ];
        
        % Check each neighbor
        for n = 1:size(neighbors, 1)
            ni = neighbors(n,1);  % neighbor row
            nj = neighbors(n,2);  % neighbor column
            
            % Check if neighbor is within bounds and is traversable (=1 in BW)
            if ni > 0 && ni <= size(BW,1) && ...
               nj > 0 && nj <= size(BW,2) && ...
               BW(ni,nj) == 1
                
                count = count + 1;
                % Convert 2D indices to 1D for sparse matrix
                I(count) = (i-1)*size(BW,2) + j;        % current node
                J(count) = (ni-1)*size(BW,2) + nj;      % neighbor node
                
                % Scale true Euclidean distance
                dx = abs(ni - i);
                dy = abs(nj - j);
                V(count) = sqrt(dx^2 + dy^2) / scaleFactor;
            end
        end
    end
    
    % Trim excess preallocated space
    I = I(1:count);
    J = J(1:count);
    V = V(1:count);
    
    % Create sparse matrix
    G = sparse(I, J, V, size(BW,1)*size(BW,2), size(BW,1)*size(BW,2));
end
