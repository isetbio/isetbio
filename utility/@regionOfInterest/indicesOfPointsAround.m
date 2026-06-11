function idx = indicesOfPointsAround(obj, points, pointsPerSmaple, samplingPoints, maxDistance)

    roiOutline = obj.outline;
    nSamples = numel(roiOutline.x);
    sampledOutline(:,1) = interp1(1:nSamples, roiOutline.x, linspace(1,nSamples, samplingPoints));
    sampledOutline(:,2) = interp1(1:nSamples, roiOutline.y, linspace(1,nSamples, samplingPoints));
    
    [D,I] = pdist2(points,sampledOutline,'euclidean','Smallest',pointsPerSmaple);
    
    idx = find(D<=maxDistance);
    idx = unique(I(idx));
end
