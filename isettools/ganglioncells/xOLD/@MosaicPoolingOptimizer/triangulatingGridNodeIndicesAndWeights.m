function [triangulatingModelRGCIndices, triangulatingModelRGCWeights, triangulatingGridNodeIndices] = ...
    triangulatingGridNodeIndicesAndWeights(theRGCposition, positionsOfModelRGCs, ...
                                           targetModelRGCindices, targetGridNodeIndices)

    distancesToModelRGCs = sqrt(sum((bsxfun(@minus, positionsOfModelRGCs, theRGCposition)).^2,2));
    [distancesToModelRGCs, idx] = sort(distancesToModelRGCs, 'ascend');
    nearestTargetModelRGCindices = targetModelRGCindices(idx);
    triangulatingGridNodeIndices = targetGridNodeIndices(idx);


    if (distancesToModelRGCs(1) == 0)
        % If the minimum distance to one of the fitted RTVF objects is
        % zero, just use that RTVF object
        triangulatingModelRGCIndices = nearestTargetModelRGCindices(1);
        triangulatingModelRGCWeights(1) = 1;
        triangulatingGridNodeIndices = triangulatingGridNodeIndices(1);
    else
        % Weights inversely proportional to the cell's distance to the
        % pnearby model RGCs
        nearbyFittedModelsNum = min([3 numel(nearestTargetModelRGCindices)]);
        triangulatingModelRGCIndices = nearestTargetModelRGCindices(1:nearbyFittedModelsNum);
        triangulatingModelRGCWeights = 1./distancesToModelRGCs(1:nearbyFittedModelsNum);
        triangulatingGridNodeIndices = triangulatingGridNodeIndices(1:nearbyFittedModelsNum);
    end
    triangulatingModelRGCWeights = triangulatingModelRGCWeights/sum(triangulatingModelRGCWeights);

end