function connectRGCsToConesBasedOnLocalDensities(obj)
% For each RGC try to connect to N cones that are not further than 1
% RGC separation away and that are not already connected to another RGC
% N is the local cone-to-RGC density ratio. The obj.coneConnectivityMatrix 
% sparse matrix is initialized here with the connections established at this
% step.

    % Compute the cone-to-RGC density ratios map at the current RGCRFpos
    densityRatiosAllRGCs = coneToRGCDensityRatiosMap(obj);

    % Indices for constructing the coneConnectivityMatrix sparse matrix
    nearestRGCindices = [];
    connectedConeIndices = [];

    activatedCones = [obj.inputConeMosaic.mConeIndices(:); obj.inputConeMosaic.lConeIndices(:)];
    [connectedConeIndices, nearestRGCindices] = doIt(obj, ...
        activatedCones,  densityRatiosAllRGCs, connectedConeIndices, nearestRGCindices);
    
    % Generate [conesNum x rgcsNum] sparse connectivity matrix
    conesNum = size(obj.inputConeMosaic.coneRFpositionsMicrons,1);
    rgcsNum = size(obj.RGCRFpositionsMicrons,1);
    weights = ones([1 numel(connectedConeIndices)]);
    obj.coneConnectivityMatrix = sparse(...
        connectedConeIndices, nearestRGCindices, weights, ...
        conesNum, rgcsNum);

    % Update centroids
    obj.updateCentroidsFromInputs(unique(nearestRGCindices));
    
end



function [connectedConeIndices, nearestRGCindices] = doIt(obj, ...
        activatedConeIndices,  densityRatiosAllRGCs, connectedConeIndices, nearestRGCindices)

    rgcsNum = size(obj.RGCRFpositionsMicrons,1);
    
    for iRGC = 1:rgcsNum

        % Find the localConeToRGCDensityRatios closest cones to each RGCRF
        [distances, idx] = RGCconnector.pdist2(...
            obj.inputConeMosaic.coneRFpositionsMicrons(activatedConeIndices,:), ...
            obj.RGCRFpositionsMicrons(iRGC,:), ...
            '', ...
            'smallest', floor(densityRatiosAllRGCs(iRGC)));

        closestConeIndices = activatedConeIndices(idx);
        if (isempty(closestConeIndices))
            continue;
        end

        % Find which of these closest cones are not already connected to
        % another RGC and also not more than 1.0 x local RGC separation
        idx = find(...
            (~ismember(closestConeIndices, connectedConeIndices)) & ...
            (distances <= obj.RGCRFspacingsMicrons(iRGC)));

        if (isempty(idx))
            continue;
        end
        closestConeIndices = closestConeIndices(idx);

        % Accumulate indices for sparse array construction 
        connectedConeIndices = cat(1, connectedConeIndices, closestConeIndices);
        nearestRGCindices = cat(1, nearestRGCindices, iRGC*ones(numel(closestConeIndices),1));

        % Update the RGC's centroid based on its inputs and weights = 1
        inputConePositions = obj.inputConeMosaic.coneRFpositionsMicrons(closestConeIndices,:);
        inputConeWeights = ones(numel(closestConeIndices),1);
       
    end % iRGC
end