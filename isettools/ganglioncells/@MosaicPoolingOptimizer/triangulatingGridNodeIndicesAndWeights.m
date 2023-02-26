function [triangulatingGridNodeIndices, triangulatingGridNodeWeights] = ...
    triangulatingGridNodeIndicesAndWeights(theRGCposition, positionsOfModelRGCs, indicesOfModelRGCs)

    distancesToModelRGCs = sqrt(sum((bsxfun(@minus, positionsOfModelRGCs, theRGCposition)).^2,2));
    [distancesToModelRGCs, idx] = sort(distancesToModelRGCs, 'ascend');
    nearestModelRGCindices = indicesOfModelRGCs(idx);

    if (distancesToModelRGCs(1) == 0)
        % If the minimum distance to one of the fitted RTVF objects is
        % zero, just use that RTVF object
        triangulatingGridNodeIndices = nearestModelRGCindices(1);
        triangulatingGridNodeWeights(1) = 1;
    else
        % Weights inversely proportional to the cell's distance to the
        % pnearby model RGCs
        nearbyFittedModelsNum = 3;
        triangulatingGridNodeIndices = nearestModelRGCindices(1:nearbyFittedModelsNum);
        triangulatingGridNodeWeights = 1./distancesToModelRGCs(1:nearbyFittedModelsNum);
    end
    triangulatingGridNodeWeights = triangulatingGridNodeWeights/sum(triangulatingGridNodeWeights);
end