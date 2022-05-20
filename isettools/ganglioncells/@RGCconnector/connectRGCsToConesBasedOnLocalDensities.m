function connectRGCsToConesBasedOnLocalDensities(obj)

    % Compute the cone-to-RGC density ratios map at the current RGCRFpos
    densityRatiosMap = coneToRGCDensityRatiosMap(obj);

    % Cell array with connected cones to each RGC
    rgcsNum = size(obj.RGCRFpositionsMicrons,1);
    obj.RGCRFinputs = cell(1, rgcsNum);
    obj.RGCRFweights = cell(1, rgcsNum);

    % List of cones already connected
    activatedCones = [obj.inputConeMosaic.mConeIndices(:); obj.inputConeMosaic.lConeIndices(:)];
    connectedConeIndices = [];
    connectedConeIndices = connectCones(obj, activatedCones, connectedConeIndices, densityRatiosMap);
    
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

function connectedConeIndices = connectCones(obj, activatedConeIndices, connectedConeIndices, densityRatiosMap)

    rgcsNum = size(obj.RGCRFpositionsMicrons,1);
    for iRGC = 1:rgcsNum
        fprintf('Working on RGC %d of %d\n', iRGC, rgcsNum);
        % Find the localConeToRGCDensityRatios closest cones to each RGCRF
        [~, idx] = RGCconnector.pdist2(...
            obj.inputConeMosaic.coneRFpositionsMicrons(activatedConeIndices,:), ...
            obj.RGCRFpositionsMicrons(iRGC,:), ...
            '', ...
            'smallest', floor(densityRatiosMap(iRGC)));

        closestConeIndices = activatedConeIndices(idx);

        % Find which of these closest cones are and not already connected
        idx = find(~ismember(closestConeIndices, connectedConeIndices));
     
        % Only connect these cones to the RGC
        obj.RGCRFinputs{iRGC} = closestConeIndices(idx);

        % Since these are the non-overlapping cones, w = 1
        obj.RGCRFweights{iRGC} = ones(1,numel(idx));

        % Keep track of connected cone indices
        connectedConeIndices = cat(1, connectedConeIndices, obj.RGCRFinputs{iRGC});
    end
end