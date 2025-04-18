clear all
close all
clc;
%%
map_rgb = imread("map.pgm");
% BW = im2bw(map_rgb,0.9);
BW1 = imbinarize(im2gray(map_rgb), 0.65);


se = strel('disk', 5);  % 5-pixel radius expansion
BW = ~imdilate(~BW1, se);

% plot_map4(BW)

map_size = size(BW);
figure; imshow(BW); hold on;
title('Select 1 GOAL point');

[u,v]=world2pixel(-3,1);
% start_pos=x_start;
start=[u,v];

% Select goal point
[goal_x, goal_y] = ginput(1); % Second click
goal = [round(goal_x), round(goal_y)];
plot(goal(1), goal(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');   % Plot goal point

%%

function [u,v] = world2pixel(x,y)
    u=20*floor(x)+200;
    v=-20*floor(y)+184;
end



function [x,y] = pixel2world(u,v)
    
    x=(u-200)/20;
    y=(v-184)/(-20);
    
end

