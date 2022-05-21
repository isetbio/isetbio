function connectRGCsToConesBasedOnLocalDensities(obj)

    % Compute the cone-to-RGC density ratios map at the current RGCRFpos
    densityRatiosMap = coneToRGCDensityRatiosMap(obj);

    % Indices for constructing the coneConnectivityMatrix sparse matrix
    nearestRGCindices = [];
    connectedConeIndices = [];

    activatedCones = [obj.inputConeMosaic.mConeIndices(:); obj.inputConeMosaic.lConeIndices(:)];
    [connectedConeIndices, nearestRGCindices] = connectCones(obj, ...
        activatedCones,  densityRatiosMap, connectedConeIndices, nearestRGCindices);
    
    % Generate [conesNum x rgcsNum] sparse connectivity matrix
    conesNum = size(obj.inputConeMosaic.coneRFpositionsMicrons,1);
    rgcsNum = size(obj.RGCRFpositionsMicrons,1);
    weights = ones([1 numel(connectedConeIndices)]);

    obj.coneConnectivityMatrix = sparse(...
        connectedConeIndices, nearestRGCindices, weights, ...
        conesNum, rgcsNum);


    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 1, ...
       'heightMargin',  0.03, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.05, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.0);

    % Visualize current connectivity
    hFig = figure(997); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1650 1000]);

    ax = subplot('Position', subplotPosVectors(1,1).v);
    [~,~,XLims, YLims] = obj.visualizeInputMosaics(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'thetaSamples', 30, ...
        'titleString', 'starting positions');
    set(ax, 'FontSize', 16)

    ax = subplot('Position', subplotPosVectors(2,1).v);
    obj.visualizeConnectivity(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'XLims', XLims, ...
        'YLims', YLims);

end

function [connectedConeIndices, nearestRGCindices] = connectCones(obj, ...
        activatedConeIndices,  densityRatiosMap, connectedConeIndices, nearestRGCindices)

    rgcsNum = size(obj.RGCRFpositionsMicrons,1);
    
    for iRGC = 1:rgcsNum
        % Find the localConeToRGCDensityRatios closest cones to each RGCRF
        [~, idx] = RGCconnector.pdist2(...
            obj.inputConeMosaic.coneRFpositionsMicrons(activatedConeIndices,:), ...
            obj.RGCRFpositionsMicrons(iRGC,:), ...
            '', ...
            'smallest', floor(densityRatiosMap(iRGC)));

        closestConeIndices = activatedConeIndices(idx);
        if (isempty(closestConeIndices))
            continue;
        end

        % Find which of these closest cones are and not already connected
        idx = find(~ismember(closestConeIndices, connectedConeIndices));
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
        [~, obj.RGCRFcentroidsFromInputs(iRGC,:)] = var(inputConePositions,inputConeWeights,1);
    end
end