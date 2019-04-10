function clout_out = ptcloud_edge_filter(cloud_in)
% extract edges from an organized RGB-D point cloud
%   Author: Maani Ghaffari Jadidi
%   Date:   January 1, 2019

BW = edge(rgb2gray(cloud_in.Color),'canny');
X = reshape(cloud_in.Location,[], 3, 1);
C = reshape(cloud_in.Color,[], 3, 1);
X = X(reshape(BW, [], 1), :);
C = C(reshape(BW, [], 1), :);
C = C(~isnan(X(:,1)),:);
X = X(~isnan(X(:,1)),:);

clout_out = pointCloud(X, 'Color', C);
