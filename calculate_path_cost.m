% total_cost = calculate_path_cost(path);

function total_cost = calculate_path_cost(path_points)
    % Calculate the total Euclidean distance cost of a path
    %
    % Inputs:
    %   path_points: N x 2 matrix where each row is [x, y] coordinates of a point
    %               or N x 1 vector of linear indices (will be converted to x,y)
    %
    % Output:
    %   total_cost: Total Euclidean distance of the path
    
    % Check if input is linear indices or coordinates
    if size(path_points, 2) == 1
        % Convert linear indices to subscripts
        [rows, cols] = ind2sub([size(BW,1), size(BW,2)], path_points);
        path_points = [cols, rows];  % Convert to [x,y] format
    end
    
    total_cost = 0;
    
    % Calculate cost between consecutive points
    for i = 1:(size(path_points, 1) - 1)
        point1 = path_points(i, :);
        point2 = path_points(i + 1, :);
        
        % Calculate Euclidean distance
        segment_cost = sqrt(sum((point2 - point1).^2));
        total_cost = total_cost + segment_cost;
    end
    
    % Print detailed information
    fprintf('Path length (number of points): %d\n', size(path_points, 1));
    fprintf('Total path cost: %.2f\n', total_cost);
    fprintf('Average segment cost: %.2f\n', total_cost/(size(path_points,1)-1));
end