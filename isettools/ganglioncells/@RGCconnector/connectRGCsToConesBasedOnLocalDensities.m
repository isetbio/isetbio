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

    % connect to L+M cones only
    consideredConeIndices = [obj.inputConeMosaic.mConeIndices(:); obj.inputConeMosaic.lConeIndices(:)];
    [connectedConeIndices, nearestRGCindices] = doIt(obj, ...
        consideredConeIndices,  densityRatiosAllRGCs, connectedConeIndices, nearestRGCindices);
    
    % Generate [conesNum x rgcsNum] sparse connectivity matrix
    conesNum = size(obj.inputConeMosaic.coneRFpositionsMicrons,1);
    rgcsNum = size(obj.RGCRFpositionsMicrons,1);
    weights = ones([1 numel(connectedConeIndices)]);
    obj.coneConnectivityMatrix = sparse(...
        connectedConeIndices, nearestRGCindices, weights, ...
        conesNum, rgcsNum);

    % Update centroids
    obj.updateCentroidsFromInputs(unique(nearestRGCindices));

    % Also set the centroids of all zero-input RGCs
    ss = squeeze(sum(obj.coneConnectivityMatrix,1));
    zeroInputRGCindices = find(ss == 0);
    obj.updateCentroidsFromInputs(zeroInputRGCindices);
end



function [connectedConeIndices, nearestRGCindices] = doIt(obj, ...
        consideredConeIndices,  densityRatiosAllRGCs, connectedConeIndices, nearestRGCindices)

    rgcsNum = size(obj.RGCRFpositionsMicrons,1);
    
    for iRGC = 1:rgcsNum

        % Find the localConeToRGCDensityRatios closest cones to each RGCRF
        [distances, idx] = RGCconnector.pdist2(...
            obj.inputConeMosaic.coneRFpositionsMicrons(consideredConeIndices,:), ...
            obj.RGCRFpositionsMicrons(iRGC,:), ...
            '', ...
            'smallest', max([1 floor(densityRatiosAllRGCs(iRGC))]));

        closestConeIndices = consideredConeIndices(idx);
        if (isempty(closestConeIndices))
            continue;
        end

        % Find which of these closest cones are not already connected to
        % another RGC and also not more than 0.6 x local RGC separation
        idx = find(...
            (~ismember(closestConeIndices, connectedConeIndices)) & ...
            (distances <= 0.6*obj.RGCRFspacingsMicrons(iRGC)));

        if (isempty(idx))
            continue;
        end
        closestConeIndices = closestConeIndices(idx);

        % Accumulate indices for sparse array construction 
        connectedConeIndices = cat(1, connectedConeIndices, closestConeIndices);
        nearestRGCindices = cat(1, nearestRGCindices, iRGC*ones(numel(closestConeIndices),1));
    end % iRGC
end