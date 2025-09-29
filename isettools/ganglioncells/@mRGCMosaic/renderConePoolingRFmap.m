function hFig = renderConePoolingRFmap(obj, theRGCindex, varargin)

    % Parse optional input
    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('clearAxes', true, @islogical);
    p.addParameter('renderConeMap', true, @islogical);
    p.addParameter('renderSubregionContour', true, @islogical)
    p.addParameter('tickSeparationArcMin', 3, @isscalar);
    p.addParameter('minConeWeightIncluded', mRGCMosaic.sensitivityAtPointOfOverlap, @isscalar);
    p.addParameter('scaleBarDegs', [], @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('doNotLabelScaleBar', false, @islogical);
    p.addParameter('forceRegenerateVisualizationCache', false, @islogical);
    p.addParameter('xSupport', [], @(x)(isempty(x)||(isnumeric(x))));
    p.addParameter('ySupport', [], @(x)(isempty(x)||(isnumeric(x))));
    p.addParameter('gridless', false, @islogical);
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('noXTicks', false, @islogical);
    p.addParameter('noYTicks', false, @islogical);
    p.parse(varargin{:});

    hFig = p.Results.figureHandle;
    theAxes = p.Results.axesHandle;
    clearAxes = p.Results.clearAxes;
    renderConeMap = p.Results.renderConeMap;
    renderSubregionContour = p.Results.renderSubregionContour;
    minConeWeightIncluded = p.Results.minConeWeightIncluded;
    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    scaleBarDegs = p.Results.scaleBarDegs;
    doNotLabelScaleBar = p.Results.doNotLabelScaleBar;
    forceRegenerateVisualizationCache = p.Results.forceRegenerateVisualizationCache;
    xSupport = p.Results.xSupport;
    ySupport = p.Results.ySupport;
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    noXTicks = p.Results.noXTicks;
    noYTicks = p.Results.noYTicks;
    gridless = p.Results.gridless;

    % Generate the axes
    if (isempty(theAxes))
        hFig = figure(); clf;
        set(hFig, 'Color', [1 1 1], 'Position', [10 10 700 500]);
        theAxes = subplot('Position', [0.05 0.06 0.91 0.91]);
    end

    % Generate visualization cache for this single cell 
    centerSubregionContourSamples = 48;
    contourGenerationMethod = 'ellipseFitToPooledConeApertureImage';

    obj.generateVisualizationCache(xSupport, ySupport, centerSubregionContourSamples, contourGenerationMethod, theRGCindex, minConeWeightIncluded, ...
        'forceRegenerateVisualizationCache', forceRegenerateVisualizationCache);

    % The cell position
    theRGCpositionDegs = obj.rgcRFpositionsDegs(theRGCindex,:);

    % The cell's RF center contour
    if (renderSubregionContour)
        theContourData = obj.visualizationCache.rfCenterContourData{1};
    else
        theContourData = [];
    end

    % Retrieve this cell's surround cone indices and weights
    if (~isempty(obj.rgcRFsurroundConeConnectivityMatrix))
        surroundConnectivityVector = full(squeeze(obj.rgcRFsurroundConeConnectivityMatrix(:, theRGCindex)));
        surroundConeIndices = find(surroundConnectivityVector > max(surroundConnectivityVector)*minConeWeightIncluded);
        surroundConeWeights = surroundConnectivityVector(surroundConeIndices);

        % The surround cones, which include the center cones
        flatTopSaturationLevel = 0.4;
        mRGCMosaic.renderInputConeMosaicSubregionPoolingMap(theAxes, ...
            obj.inputConeMosaic, ...
            theContourData, ...
            surroundConeIndices, ...
            surroundConeWeights, ...
            flatTopSaturationLevel, ...
            theRGCpositionDegs, ...
            tickSeparationArcMin, ...
            scaleBarDegs, ...
            gridless, noXLabel, noYLabel, ...
            'clearAxes', clearAxes, ...
            'doNotLabelScaleBar', doNotLabelScaleBar, ...
            'renderConeMap', renderConeMap, ...
            'plotTitle', sprintf('min cone weight visualized: %2.3f', minConeWeightIncluded));

    else
        connectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, theRGCindex)));
        centerConeIndices = find(connectivityVector > minConeWeightIncluded);
        centerConeWeights = connectivityVector(centerConeIndices);

        % The center cones
        flatTopSaturationLevel = 1.0;
        mRGCMosaic.renderInputConeMosaicSubregionPoolingMap(theAxes, ...
            obj.inputConeMosaic, ...
            theContourData, ...
            centerConeIndices, ...
            centerConeWeights, ...
            flatTopSaturationLevel, ...
            theRGCpositionDegs, ...
            tickSeparationArcMin, ...
            scaleBarDegs, ...
            gridless, noXLabel, noYLabel, ...
            'clearAxes', clearAxes, ...
            'doNotLabelScaleBar', doNotLabelScaleBar, ...
            'renderConeMap', renderConeMap, ...
            'plotTitle', sprintf('min cone weight visualized: %2.3f', minConeWeightIncluded));
    end
end