function ptCloudOut = pcRangeFilter(ptCloudIn, max_range, min_range)
% Remove all the points in ptCloudIn with range larger than threshold

% Calculate range for all points
pt_location = ptCloudIn.Location;
pt_color = ptCloudIn.Color;
pt_range = sqrt(sum(pt_location .* pt_location, 2));

% Filtering
pt_location(pt_range > max_range | pt_range < min_range, :) = [];
pt_color(pt_range > max_range | pt_range < min_range, :) = [];
ptCloudOut = pointCloud(pt_location, 'Color', pt_color);

end
