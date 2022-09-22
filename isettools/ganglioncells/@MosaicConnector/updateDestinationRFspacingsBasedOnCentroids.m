function updateDestinationRFspacingsBasedOnCentroids(obj)

    % Find all the distances
    [distancesToNearbyDestinationRFs, nearbyDestinationRFIndices] = ...
        MosaicConnector.pdist2(...
            obj.destinationRFcentroidsFromInputs, ...
            obj.destinationRFcentroidsFromInputs, ...
            'smallest', obj.wiringParams.maxNeighborsNum+1);

    % Initialize destination RF spacings from their current centroids 
    obj.destinationRFspacingsFromCentroids = inf(1,size(obj.destinationRFcentroidsFromInputs,1));
    
    nonZeroInputDestinationRFindices = find(~isinf(obj.destinationRFcentroidsFromInputs(:,1)));
    spacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(obj.destinationRFcentroidsFromInputs);

    maxNeighborNormDistance = obj.wiringParams.maxNeighborNormDistance;
    destinationRFspacingsFromCentroids = zeros(1,numel(nonZeroInputDestinationRFindices));

    % Update the spacings for all destinationRFs with at least 1 input
    parfor iRF = 1:numel(nonZeroInputDestinationRFindices)
         theDestinationRFindex = nonZeroInputDestinationRFindices(iRF);

         nearbyDestinationRFindicesForThisDestinationRF = nearbyDestinationRFIndices(:,theDestinationRFindex);
         distancesToNearbyDestinationRFsForThisDestinationRF = distancesToNearbyDestinationRFs(:,theDestinationRFindex);

         % Exclude nearbyRGCs that are further than a maxDistance
         maxDistance = maxNeighborNormDistance * spacingsMicrons(theDestinationRFindex);

         nearbyDestinationRFindicesForThisDestinationRF = ...
             nearbyDestinationRFindicesForThisDestinationRF(find(distancesToNearbyDestinationRFsForThisDestinationRF < maxDistance));

         destinationRFspacingsFromCentroids(iRF) = mean(spacingsMicrons(nearbyDestinationRFindicesForThisDestinationRF));
    end

    obj.destinationRFspacingsFromCentroids(nonZeroInputDestinationRFindices) = destinationRFspacingsFromCentroids;
end
