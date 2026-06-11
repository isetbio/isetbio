function idx = indicesOfPointsInside(obj, points)

    assert(size(points,2) == 2, 'points must be an [Nx2] matrix');
     
    % Compute the roiOutline
    roiOutline = obj.outline;
    
    % Find indices of cones within the ROI outline
    [in,on] = inpolygon( points(:,1), points(:,2),...
                         roiOutline.x, roiOutline.y);
    idx = find((in|on));
end
