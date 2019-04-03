%   Author: Maani Ghaffari Jadidi
%   Date:   January 1, 2019

clc; clear; close all

% Register Two Point Clouds
dataFile = fullfile(toolboxdir('vision'), 'visiondata', 'livingRoom.mat');
load(dataFile);

ptCloudRef = ptcloud_edge_filter(livingRoomData{1});
ptCloudCurrent = ptcloud_edge_filter(livingRoomData{2});


gridSize = 0.1;
fixed = pcdownsample(ptCloudRef, 'gridAverage', gridSize);
moving = pcdownsample(ptCloudCurrent, 'gridAverage', gridSize);

% make rkhs registration object
rkhs_se3 = rkhs_se3_registration();
rkhs_se3.set_ptclouds(fixed, moving);
tic
rkhs_se3.align();
toc
fprintf('\n')
tform = rkhs_se3.tform;

% Transform the current point cloud to the reference coordinate system
% defined by the first point cloud.
ptCloudAligned = pctransform(livingRoomData{2},tform);

mergeSize = 0.015;
ptCloudScene = pcmerge(livingRoomData{1}, ptCloudAligned, mergeSize);

% Stitch a Sequence of Point Clouds
% To compose a larger 3-D scene, repeat the same procedure as above to
% process a sequence of point clouds. Use the first point cloud to
% establish the reference coordinate system. Transform each point cloud to
% the reference coordinate system. This transformation is a multiplication
% of pairwise transformations.

% Store the transformation object that accumulates the transformation.
accumTform = tform; 

figure
hAxes = pcshow(ptCloudScene, 'VerticalAxis','Y', 'VerticalAxisDir', 'Down');
title('Updated world scene')
% Set the axes property for faster rendering
hAxes.CameraViewAngleMode = 'auto';
hScatter = hAxes.Children;

for i = 3:length(livingRoomData)
    ptCloudCurrent = ptcloud_edge_filter(livingRoomData{i});
       
    % Use previous moving point cloud as reference.
    fixed = moving;
    moving = pcdownsample(ptCloudCurrent, 'gridAverage', gridSize);
    
    % Apply RKHS registration.
    rkhs_se3.set_ptclouds(fixed, moving); 
    tic
    rkhs_se3.align();
    toc
    fprintf('\n')
    tform = rkhs_se3.tform;

    % Transform the current point cloud to the reference coordinate system
    % defined by the first point cloud.
    accumTform = affine3d(tform.T * accumTform.T);
    ptCloudAligned = pctransform(livingRoomData{i}, accumTform);
    
    % Update the world scene.
    ptCloudScene = pcmerge(ptCloudScene, ptCloudAligned, mergeSize);

    % Visualize the world scene.
    hScatter.XData = ptCloudScene.Location(:,1);
    hScatter.YData = ptCloudScene.Location(:,2);
    hScatter.ZData = ptCloudScene.Location(:,3);
    hScatter.CData = ptCloudScene.Color;
    drawnow('limitrate')
end

% During the recording, the Kinect was pointing downward. To visualize the
% result more easily, let's transform the data so that the ground plane is
% parallel to the X-Z plane.
angle = -pi/10;
A = [1,0,0,0;...
     0, cos(angle), sin(angle), 0; ...
     0, -sin(angle), cos(angle), 0; ...
     0 0 0 1];
ptCloudScene = pctransform(ptCloudScene, affine3d(A));
pcshow(ptCloudScene, 'VerticalAxis','Y', 'VerticalAxisDir', 'Down', ...
        'Parent', hAxes)
title('Updated world scene')
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
