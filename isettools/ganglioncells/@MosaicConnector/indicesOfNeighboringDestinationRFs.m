function [nearbyDestinationRFIndices,distancesToNearbyDestinationRFs] = indicesOfNeighboringDestinationRFs(obj, destinationRFIndex, varargin)
       
    p = inputParser;
    p.addParameter('maxNeighborsNum', [], @(x)((isempty(x))||(isscalar(x))));
    p.addParameter('ignoreInfCentroids', false, @islogical);
    p.parse(varargin{:});
    
    maxNormDistance = obj.wiringParams.maxNeighborNormDistance;
    if (isempty(p.Results.maxNeighborsNum))
        maxNeighborsNum = obj.wiringParams.maxNeighborsNum;
    else
        maxNeighborsNum = p.Results.maxNeighborsNum;
    end

    % Find the indices of the neighboring destination RFs
    [distancesToNearbyDestinationRFs, nearbyDestinationRFIndices] = MosaicConnector.pdist2(...
            obj.destinationRFcentroidsFromInputs, ...
            obj.destinationRFcentroidsFromInputs(destinationRFIndex,:), ...
            'smallest', maxNeighborsNum+1);
                    
    % Exclude the destinationRFIndex, which has 0 distance
    [nearbyDestinationRFIndices,ia] = setdiff(nearbyDestinationRFIndices, destinationRFIndex);
    distancesToNearbyDestinationRFs = distancesToNearbyDestinationRFs(ia);

    % Exclude nearbyRGCs that are further than a maxDistance
    if (~p.Results.ignoreInfCentroids) && (isinf(obj.destinationRFcentroidsFromInputs(destinationRFIndex)))
        error('local spacing should not be inf here')
    end
    maxDistance = maxNormDistance * obj.destinationRFspacingsFromCentroids(destinationRFIndex);
    
    idx = find(distancesToNearbyDestinationRFs <= maxDistance);
    nearbyDestinationRFIndices = nearbyDestinationRFIndices(idx);
    distancesToNearbyDestinationRFs = distancesToNearbyDestinationRFs(idx);
    
    % Sort them according to their distance to the destinationRFIndex
    [~,idx] = sort(distancesToNearbyDestinationRFs, 'ascend');
    nearbyDestinationRFIndices = nearbyDestinationRFIndices(idx);
    distancesToNearbyDestinationRFs = distancesToNearbyDestinationRFs(idx);
end