function obj = cropToSizeAtEccentricity(obj, sizeDegs, eccentricityDegs, varargin)
    % Parse input
    p = inputParser;
    p.addRequired('sizeDegs', @isnumeric);
    p.addRequired('eccentricityDegs', @isnumeric);
    p.addParameter('name', '', @ischar);
    p.addParameter('visualizeSpatialRelationshipToSourceMosaic', false, @islogical);
    p.addParameter('extraSupportDegsForInputConeMosaic', 0, @isscalar);
    p.addParameter('beVerbose', false, @islogical);

    % Parse input
    p.parse(sizeDegs, eccentricityDegs, varargin{:});
    assert(isempty(eccentricityDegs)||numel(eccentricityDegs)==2, 'eccentricity must be either [] OR a 2-element vector');
    assert(isempty(sizeDegs)||numel(sizeDegs)==2, 'size must be either [] (no cropping) OR a 2-element vector');

    if (isempty(sizeDegs))
        % No cropping
        return;
    end
    
    if (isempty(eccentricityDegs))
        % Crop at the mosaic's center
        eccentricityDegs = obj.eccentricityDegs;
    end

    name = p.Results.name;
    visualizeSpatialRelationshipToSourceMosaic = p.Results.visualizeSpatialRelationshipToSourceMosaic;
    beVerbose = p.Results.beVerbose;
    extraSupportDegsForInputConeMosaic = p.Results.extraSupportDegsForInputConeMosaic;

    % Define theROI based on the passed eccentricity and size
    theROI = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', eccentricityDegs, ...
            'width', sizeDegs(1), ...
            'height', sizeDegs(2), ...
            'rotation', 0.0...
        ));

    if (visualizeSpatialRelationshipToSourceMosaic)
        [hFig, ax] = obj.visualize();
        hold(ax, 'on');
        theROI.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax);
        set(ax, 'XLim', obj.eccentricityDegs(1) + 0.5*obj.sizeDegs(1)*[-1 1], ...
                'YLim', obj.eccentricityDegs(2) + 0.5*obj.sizeDegs(2)*[-1 1], ...
                'XTick', eccentricityDegs(1), ...
                'YTick', eccentricityDegs(2), ...
                'XTickLabel', sprintf('%2.1f', eccentricityDegs(1)), ...
                'YTickLabel', sprintf('%2.1f', eccentricityDegs(2)));
    end


    % Find the indices of the RGCs whose position is within theROI
    keptRGCindices = theROI.indicesOfPointsInside(obj.rgcRFpositionsDegs);
    if (numel(keptRGCindices) == 0)
        error('No RGCs left after cropping to (%2.1f x %2.1f) at (%2.1f,%2.1f) degs\n', ...
            sizeDegs(1), sizeDegs(2), eccentricityDegs(1), eccentricityDegs(2));
        return;
    end

    % Find cone indices to keep based on the surrounds of the kept RGCs
    keptConeIndices = [];
    if (~isempty(obj.rgcRFcenterConeConnectivityMatrix))
        for iRGC = 1:numel(keptRGCindices)
            coneIndices = find(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, keptRGCindices(iRGC)))>0.000001);
            keptConeIndices = unique(cat(1,keptConeIndices(:),coneIndices(:)));
        end

        if (extraSupportDegsForInputConeMosaic > 0)
            minConePos = min(obj.inputConeMosaic.coneRFpositionsDegs(keptConeIndices,:), [], 1);
            maxConePos = max(obj.inputConeMosaic.coneRFpositionsDegs(keptConeIndices,:), [], 1);
            keptConeIndices = find(...
                (obj.inputConeMosaic.coneRFpositionsDegs(:,1) >= minConePos(1)-extraSupportDegsForInputConeMosaic) & ...
                (obj.inputConeMosaic.coneRFpositionsDegs(:,1) <= maxConePos(1)+extraSupportDegsForInputConeMosaic) & ...
                (obj.inputConeMosaic.coneRFpositionsDegs(:,2) >= minConePos(2)-extraSupportDegsForInputConeMosaic) & ...
                (obj.inputConeMosaic.coneRFpositionsDegs(:,2) <= maxConePos(2)+extraSupportDegsForInputConeMosaic));
        end
    else
        error('@mRGCMosaic does not contain connectivity matrices');
    end

    % Determine the indices of cones to keep
    xPos = obj.inputConeMosaic.coneRFpositionsDegs(keptConeIndices,1);
    yPos = obj.inputConeMosaic.coneRFpositionsDegs(keptConeIndices,2);
    xMinPos = min(xPos(:));
    xMaxPos = max(xPos(:));
    yMinPos = min(yPos(:));
    yMaxPos = max(yPos(:));

    keptAllConeTypeIndices = find(...
        (obj.inputConeMosaic.coneRFpositionsDegs(:,1) >= xMinPos) & ...
        (obj.inputConeMosaic.coneRFpositionsDegs(:,1) <= xMaxPos) & ...
        (obj.inputConeMosaic.coneRFpositionsDegs(:,2) >= yMinPos) & ...
        (obj.inputConeMosaic.coneRFpositionsDegs(:,2) <= yMaxPos) ...
        );

    % Update obj.inputConeMosaic to only include the keptConeIndices
    obj.inputConeMosaic.cropMosaicToIncluceConesWithIndices(keptAllConeTypeIndices);

    % Update the rgcRFpositions and spacings of the cropped mRGCMosaic 
    obj.rgcRFpositionsMicrons = obj.rgcRFpositionsMicrons(keptRGCindices,:);
    obj.rgcRFpositionsDegs = obj.rgcRFpositionsDegs(keptRGCindices,:);
    obj.rgcRFspacingsMicrons = obj.rgcRFspacingsMicrons(keptRGCindices);
    obj.rgcRFspacingsDegs = obj.rgcRFspacingsDegs(keptRGCindices);
    obj.rgcsNum = numel(keptRGCindices);

    % Crop obj.rgcRFcenterConeConnectivityMatrix to generate croppedCenterConePoolingMatrix
    obj.rgcRFcenterConeConnectivityMatrix = obj.rgcRFcenterConeConnectivityMatrix(keptAllConeTypeIndices, keptRGCindices);

    if (~isempty(obj.rgcRFsurroundConeConnectivityMatrix))
        % Crop obj.rgcRFsurroundConeConnectivityMatrix
        obj.rgcRFsurroundConeConnectivityMatrix = obj.rgcRFsurroundConeConnectivityMatrix(keptAllConeTypeIndices, keptRGCindices);
    end

    if (~isempty(obj.rgcRFcenterConeConnectivityMatrixNoOverlap))
        % Crop obj.rgcRFsurroundConeConnectivityMatrixNoOverlap
        obj.rgcRFcenterConeConnectivityMatrixNoOverlap = obj.rgcRFcenterConeConnectivityMatrixNoOverlap(keptAllConeTypeIndices, keptRGCindices);
    end

    if (~isempty(obj.exclusivelyConnectedInputConeIndicesNum))
        obj.exclusivelyConnectedInputConeIndicesNum = obj.exclusivelyConnectedInputConeIndicesNum(keptRGCindices);
    end

    if (~isempty(obj.responseGains))
        obj.responseGains = obj.responseGains(keptRGCindices);
    end

    % Update the eccentricity and size of the cropped mRGCMosaic 
    minRFpositionDegs = squeeze(min(obj.rgcRFpositionsDegs,[],1));
    maxRFpositionDegs = squeeze(max(obj.rgcRFpositionsDegs,[],1));
    minRFpositionMicrons = squeeze(min(obj.rgcRFpositionsMicrons,[],1));
    maxRFpositionMicrons = squeeze(max(obj.rgcRFpositionsMicrons,[],1));

    obj.eccentricityDegs = 0.5*(maxRFpositionDegs+minRFpositionDegs);
    obj.eccentricityMicrons = 0.5*(maxRFpositionMicrons+minRFpositionMicrons);
 
    obj.sizeDegs = maxRFpositionDegs-minRFpositionDegs;

    % Reset the visualizationCache
    obj.visualizationCache = [];

    if (visualizeSpatialRelationshipToSourceMosaic)
        hFig = figure(55);
        set(hFig, 'Position', [10 10 1700 800]);
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
            'rowsNum', 1, ...
            'colsNum', 2, ...
            'heightMargin',  0.06, ...
            'widthMargin',    0.05, ...
            'leftMargin',     0.04, ...
            'rightMargin',    0.00, ...
            'bottomMargin',   0.04, ...
            'topMargin',      0.0);

        ax = subplot('position', subplotPosVectors(1,1).v);
        obj.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax);

        ax = subplot('position', subplotPosVectors(1,2).v);
        obj.inputConeMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax);
    end

end
