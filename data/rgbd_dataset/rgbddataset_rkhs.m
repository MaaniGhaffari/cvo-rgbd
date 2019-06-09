function rgbddataset_rkhs(dataset_name, varargin)


file_folder = ...
    strcat(dataset_name, '/pcd_ds/');
file_list = dir(strcat(file_folder, '*.pcd'));
file_number = length(file_list);

% parameters
option = [];
option.max_range = 4;
option.min_range = 0.8;

if nargin > 1 && ~isempty(varargin{1})
    option.gridSize = varargin{1};
else
    option.gridSize = 0.05;
end

if nargin > 2 && strcmp(varargin{2}, 'quiet')
    option.quiet = true;
else
    option.quiet = false;
end

if nargin > 3 && ~isempty(varargin{3})
    option.numiter = varargin{3};
else
    option.numiter = file_number;
end

result = cell(file_number,1);

target = pcread(file_list(1).name);
% ptCloudRef = ptcloud_edge_filter(target);
ptCloudRef = pcRangeFilter(target, option.max_range, option.min_range);

% downsample point clouds using grids
fixed = pcdownsample(ptCloudRef, 'gridAverage', option.gridSize);

tform = affine3d;
result{1} = tform;
registration_time = zeros(file_number-1,1);

% make rkhs registration object
rkhs_se3 = rkhs_se3_registration();

for i = 2:option.numiter
    try
        source = pcread(file_list(i).name);
%         ptCloudCurrent = ptcloud_edge_filter(source);
        ptCloudCurrent = pcRangeFilter(source, option.max_range, option.min_range);

        % downsample point clouds using grids
        moving = pcdownsample(ptCloudCurrent, 'gridAverage', option.gridSize);
        
        if ~option.quiet
            disp(['Registering: ', file_list(i-1).name, ' and ', file_list(i).name])
            disp(['Registering: ', num2str(i-1), ' and ', num2str(i)])
        end
        
        t0 = tic;
        rkhs_se3.set_ptclouds(fixed, moving);
        rkhs_se3.align();
        tform = rkhs_se3.tform;
        registration_time(i-1) = toc(t0);
        if ~option.quiet
            disp(['Registration time in seconds: ', num2str(registration_time(i-1))])
        end
        
        fixed = moving;
        result{i} = tform;
    catch
        if ~option.quiet
            disp(['Registration failed at between frames: ', num2str(i-1), ' and ', num2str(i)])
        end
        
        registration_time(i-1) = nan;
        fixed = moving;
        result{i} = nan;
    end
    if ~option.quiet
        disp('------------')
    end
end

file_name = strcat(dataset_name, '_', datestr(now, 'dd-mmm-yyyy-HH-MM-SS'), '.mat');
save(file_name, 'result', 'option', 'dataset_name', 'registration_time');

end
