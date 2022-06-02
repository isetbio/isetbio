function updateLocalRGCRFspacingsBasedOnCurrentCentroids(obj)
    
    % Also compute local RGCRFspacings by averaging the
    % RGCRFspacingsMicrons over its neighbors
    
    maxNeighborsNum = obj.wiringParams.maxNeighborsNum;
    maxNormDistance = obj.wiringParams.maxNeighborNormDistance;

    [distancesToNearbyRGCs, nearbyRGCindices] = RGCconnector.pdist2(...
            obj.RGCRFcentroidsFromInputs, ...
            obj.RGCRFcentroidsFromInputs, '', ...
            'smallest', maxNeighborsNum+1);

    % Compute RGCRF spacings from their current centroids 
    obj.localRGCRFspacingsMicrons = inf(1,size(obj.RGCRFcentroidsFromInputs,1));
    
    nonZeroInputRGCIDs = find(~isinf(obj.RGCRFcentroidsFromInputs(:,1)));
    localRGCRFspacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(obj.RGCRFcentroidsFromInputs);

    for iiiRGC = 1:numel(nonZeroInputRGCIDs)
         iRGC = nonZeroInputRGCIDs(iiiRGC);
         nearbyRGCindicesForThisRGC = nearbyRGCindices(:,iRGC);
         distancesToNearbyRGCsForThisRGC = distancesToNearbyRGCs(:,iRGC);
         % Exclude nearbyRGCs that are further than a maxDistance
         maxDistance = maxNormDistance * localRGCRFspacingsMicrons(iRGC);
         nearbyRGCindicesForThisRGC = nearbyRGCindicesForThisRGC(find(distancesToNearbyRGCsForThisRGC < maxDistance));
         obj.localRGCRFspacingsMicrons(iRGC) = mean(localRGCRFspacingsMicrons(nearbyRGCindicesForThisRGC));
    end
end