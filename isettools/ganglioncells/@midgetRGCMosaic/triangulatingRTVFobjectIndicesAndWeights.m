function [triangulatingRTVFobjIndices, triangulatingRTVFobjWeights] = triangulatingRTVFobjectIndicesAndWeights(obj,iRGC)

    % Retrieve this cell's # of center cone indices
    connectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, iRGC)));
    indicesOfCenterCones = find(connectivityVector > 0.0001);

    % Retrieve the indices of the fitted RTVF objects that have the
    % same # of center cones
    centerConeMatchObjIndices = find(obj.theConesNumPooledByTheRFcenterGrid == numel(indicesOfCenterCones));

    % Compute distance based weights for this RGC and the fitted RTVF objects
    distancesToSamplingGridPositions = sqrt(sum((bsxfun(@minus, obj.theSamplingPositionGrid(centerConeMatchObjIndices,:), obj.rgcRFpositionsDegs(iRGC,:))).^2,2));
    [distancesToSamplingGridPositions, neighboringSamplingGridPositions] = sort(distancesToSamplingGridPositions, 'ascend');

    if (distancesToSamplingGridPositions(1) == 0)
            % If the minimum distance to one of the fitted RTVF objects is
            % zero, just use that RTVF object
            neighboringSamplingGridPositions = neighboringSamplingGridPositions(1);
            distancesToSamplingGridPositions = distancesToSamplingGridPositions(1);
            triangulatingRTVFobjIndices = centerConeMatchObjIndices(neighboringSamplingGridPositions);
            triangulatingRTVFobjWeights(1) = 1;
    else
            % We will sum the weigted surrounds from the nearest nearbyFittedModelsNum
            % with weights inversely proportional to the cell's distance to the position of the nearbyFittedModels
            nearbyFittedModelsNum = 3;
            neighboringSamplingGridPositions = neighboringSamplingGridPositions(1:nearbyFittedModelsNum);
            distancesToSamplingGridPositions = distancesToSamplingGridPositions(1:nearbyFittedModelsNum);
            triangulatingRTVFobjIndices = centerConeMatchObjIndices(neighboringSamplingGridPositions);
            triangulatingRTVFobjWeights = 1./distancesToSamplingGridPositions;
            triangulatingRTVFobjWeights = triangulatingRTVFobjWeights/sum(triangulatingRTVFobjWeights);
    end

end
