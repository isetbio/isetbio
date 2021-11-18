function idx = indicesOfPointsAround(obj, points, maxDistance)

    roiOutline = obj.outline;
    nSamples = numel(roiOutline.x);
    sampledOutline(:,1) = interp1(1:nSamples, roiOutline.x, linspace(1,nSamples, 100));
    sampledOutline(:,2) = interp1(1:nSamples, roiOutline.y, linspace(1,nSamples, 100));
    
    [D,I] = pdist2(points,sampledOutline,'euclidean','Smallest',1);
    
    idx = find(D<=maxDistance);
    idx = I(idx);
end
