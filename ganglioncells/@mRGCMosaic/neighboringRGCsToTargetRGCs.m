function [nearbyRGCindices, distanceMatrix] = ...
	neighboringRGCsToTargetRGCs(obj, theTargetRGCindices, maxNeighborsNum)

	% nearbyRGCindices: [maxNeighborsNum x numel(theTargetRGCindices)] matrix of indices of nearby RGCs to theTargetRGCindices
    % distanceMatrix: [maxNeighborsNum x numel(theTargetRGCindices)] matrix of distances of nearby RGCs to theTargetRGCindices

    [distanceMatrix, nearbyRGCindices] = MosaicConnector.pdist2(...
            obj.rgcRFpositionsDegs, ...
    		obj.rgcRFpositionsDegs(theTargetRGCindices,:), ...
            'smallest', maxNeighborsNum+1);
           
    distanceMatrix = distanceMatrix(2:end,:);
    nearbyRGCindices = nearbyRGCindices(2:end,:);
end