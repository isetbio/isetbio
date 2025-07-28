function [nearbyRGCindices, distancesToNearbyRGCs] = indicesOfNearbyRGCs(obj, targetRGC, varargin)
       
    p = inputParser;
    p.addParameter('maxNeighborsNum', 6, @(x)((isempty(x))||(isscalar(x))));
    p.addParameter('maxNormDistance', inf, @isscalar);
    p.parse(varargin{:});
    
    maxNeighborsNum = p.Results.maxNeighborsNum;
    maxNormDistance = p.Results.maxNormDistance;

    % Find the indices of the neighboring destination RFs
    [distancesToNearbyRGCs, nearbyRGCindices] = MosaicConnector.pdist2(...
            obj.rgcRFpositionsDegs, ...
            obj.rgcRFpositionsDegs(targetRGC,:), ...
            'smallest', maxNeighborsNum+1);
                    
    % Exclude the targetRGC, which has 0 distance
    [nearbyRGCindices,ia] = setdiff(nearbyRGCindices, targetRGC);
    distancesToNearbyRGCs = distancesToNearbyRGCs(ia);

    if isinf(maxNormDistance)
        return;
    end

    maxDistance = maxNormDistance * obj.rgcRFspacingsDegs(targetRGC);
    idx = find(distancesToNearbyRGCs  <= maxDistance);
    nearbyRGCindices = nearbyRGCindices(idx);
    distancesToNearbyRGCs = distancesToNearbyRGCs(idx);
    
    % Sort them according to their distance to the destinationRFIndex
    [~,idx] = sort(distancesToNearbyRGCs, 'ascend');
    nearbyRGCindices = nearbyRGCindices(idx);
    distancesToNearbyRGCs = distancesToNearbyRGCs(idx);
end