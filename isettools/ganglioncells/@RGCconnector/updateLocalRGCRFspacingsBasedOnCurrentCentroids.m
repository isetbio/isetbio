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
    localRGCRFspacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(obj.RGCRFcentroidsFromInputs);

    idx = find(isinf(localRGCRFspacingsMicrons));
    if (~isempty(idx))
        error('Inf local spacing found')
        localRGCRFspacingsMicrons(idx) = obj.RGCRFspacingsMicrons(idx);
    end

    for iRGC = 1:numel(obj.RGCRFspacingsMicrons)
         nearbyRGCindicesForThisRGC = nearbyRGCindices(:,iRGC);
         distancesToNearbyRGCsForThisRGC = distancesToNearbyRGCs(:,iRGC);
         % Exclude nearbyRGCs that are further than a maxDistance
         maxDistance = maxNormDistance * localRGCRFspacingsMicrons(iRGC);
         nearbyRGCindicesForThisRGC = nearbyRGCindicesForThisRGC(find(distancesToNearbyRGCsForThisRGC <= maxDistance));
         obj.localRGCRFspacingsMicrons(iRGC) = mean(localRGCRFspacingsMicrons(nearbyRGCindicesForThisRGC));
     end


end