function rgbddataset_cdf_plots(dataset_name)

switch dataset_name
    case 'freiburg1_desk'
        load('freiburg1_desk-gt-pose.mat');
        load('freiburg1_desk_07-May-2019-02-35-00.mat');
        load('opencv_freiburg1_desk.mat');
        %
    case 'freiburg1_desk2'
        %
    case 'freiburg1_room'
        %
    case 'freiburg2_pioneer_360'
        %
    case 'freiburg2_pioneer_slam'
        %
    case 'freiburg2_desk'
        %
    case 'freiburg3_nostructure_notexture_far'
        %
    case 'freiburg3_nostructure_texture_far'
        %
    case 'freiburg3_structure_notexture_far'
        %
    case 'freiburg3_structure_notexture_near'
        %
    case 'freiburg3_structure_texture_far'
        %
    case 'freiburg3_structure_texture_near'
        %
    case 'freiburg3_long_office_household'
        %
    otherwise
        assert(false,'The sequence name is not available.')
end


% get 4x4 tf from ROS position and quaterion (x,y,z,w). MATLAB uses
% w,x,y,z order.
ros2pose = @(t,q) ([quat2rotm(q), t; 0 0 0 1]);
tfinv = @(T) ([T(1:3,1:3)', -T(1:3,1:3)' * T(1:3,4); 0 0 0 1]);

% find the first matching pose
k = 1;
while (assoc(1) - gt_pose(k,1)) > 0 && k < size(gt_pose,1)
    k = k + 1;
end

% so(3) and R^3 distance metric
dso3 = @(R1,R2) (norm(logm(R1 * R2'), 'fro'));
dR3 = @(R1,R2,t1,t2) (norm(t1 - R1 * R2' * t2));

rot = 0;        % rotation error
tran = 0;       % translation error
rot_cv = 0;     % rotation error
tran_cv = 0;    % translation error

t = gt_pose(k,2:4)';               % position
q = gt_pose(k,[8,5,6,7]);          % orientation (quaternion)
% get 4x4 tf
T0 = ros2pose(t, q);

j = 1; % pose estimates iterator
for i = 2:length(assoc)
    % find the next matching pose
    while (assoc(i) - gt_pose(k,1)) > 0 && k < size(gt_pose,1)
        k = k + 1;
    end
        
    t = gt_pose(k,2:4)';               % position
    q = gt_pose(k,[8,5,6,7]);          % orientation (quaternion)
    % get 4x4 tf 
    T1 = ros2pose(t, q);
    Trel = tfinv(T0) * T1;
    T0 = T1;
    
    % compute rkhs registration error
    if isa(result{i}, 'affine3d')
        H = result{i}.T';
        rot(j,1) = dso3(H(1:3,1:3), Trel(1:3,1:3));
        tran(j,1) = dR3(H(1:3,1:3), Trel(1:3,1:3), H(1:3,4), Trel(1:3,4));
    else
        rot(j,1) = nan;
        tran(j,1) = nan;
    end
    %         H = [reshape(result(i-1,6:end),3,3)', result(i-1,3:5)'; 0 0 0 1];
    
    % compute opencv benchmark error; some frames failed so we fill in nan.
    % TODO: not fair actually! cvo might have larger errors while it didn't fail.
    H = [reshape(cvrgbdposes(i-1,6:end),3,3)', cvrgbdposes(i-1,3:5)'; 0 0 0 1];
    H = tfinv(H);
    if ~all(all(H == eye(4)))
        rot_cv(j,1) = dso3(H(1:3,1:3), Trel(1:3,1:3));
        tran_cv(j,1) = dR3(H(1:3,1:3), Trel(1:3,1:3), H(1:3,4), Trel(1:3,4));
    else
        rot_cv(j,1) = nan;
        tran_cv(j,1) = nan;
    end
    j = j + 1;
end

% CDF plots
fsize = 22;
fig = figure; hold on; set(gca,'TickLabelInterpreter','latex', 'fontsize', fsize); 
h{1} = cdfplot(tran); 
set(h{1}, 'linewidth', 2.5, 'linestyle', '-')
h{2} = cdfplot(tran_cv);
set(h{2}, 'linewidth', 2.5, 'linestyle', '--')
legend([h{1}, h{2}], 'RGB-D CVO', 'OpenCV RGB-D VO', 'location', 'best','Interpreter','latex')
xlabel('Position error (m)','Interpreter','latex')
ylabel('Fraction of data','Interpreter','latex')
% title(dataset_name, 'FontWeight', 'normal', 'fontsize', fsize), grid on, axis tight
title(''), grid on, axis tight
figuresize(21,12,'cm')
image_name = ['tum_cdf_',dataset_name,'_position_error'];
print(fig, image_name, '-dpng', '-r300')

fig = figure; hold on; set(gca,'TickLabelInterpreter','latex', 'fontsize', fsize); 
h{1} = cdfplot(rot * 180/pi); 
set(h{1}, 'linewidth', 2.5, 'linestyle', '-')
h{2} = cdfplot(rot_cv * 180/pi);
set(h{2}, 'linewidth', 2.5, 'linestyle', '--')
legend([h{1}, h{2}], 'RGB-D CVO', 'OpenCV RGB-D VO', 'location', 'best','Interpreter','latex')
xlabel('Orientation error (deg)','Interpreter','latex')
ylabel('Fraction of data','Interpreter','latex')
title(''), grid on, axis tight
figuresize(21,12,'cm')
image_name = ['tum_cdf_',dataset_name,'_orientation_error'];
print(fig, image_name, '-dpng', '-r300')

end
