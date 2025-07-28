function visualizeRetinalConePoolingRFmapOfRGCwithIndex(obj, theRGCindex, varargin)
    % Parse optional input
    p = inputParser;
    p.addParameter('surroundSaturationLevel', 1, @isscalar);
    p.addParameter('withFigureFormat', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('tickSeparationArcMin', 6, @isscalar);
    p.addParameter('normalizedPeakSurroundSensitivity', 0.4, @isscalar);
    p.addParameter('theAxes', [], @(x)(isempty(x)||(iscell(x))));
    p.addParameter('reverseXDir', false, @islogical);
    p.addParameter('gridlessLineWeightingFunctions', false, @islogical);
    p.addParameter('noYLabelLineWeightingFunctions', true, @islogical);
    p.addParameter('regenerateVisualizationCache', true, @islogical);

    p.parse(varargin{:});
    surroundSaturationLevel = p.Results.surroundSaturationLevel;
    ff = p.Results.withFigureFormat;
    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    normalizedPeakSurroundSensitivity = p.Results.normalizedPeakSurroundSensitivity;
    theAxes = p.Results.theAxes;
    reverseXDir = p.Results.reverseXDir;
    gridlessLineWeightingFunctions = p.Results.gridlessLineWeightingFunctions;
    noYLabelLineWeightingFunctions = p.Results.noYLabelLineWeightingFunctions;
    regenerateVisualizationCache = p.Results.regenerateVisualizationCache;

    if (regenerateVisualizationCache)
        % Generate the visualization cache
        xSupport = [];
        ySupport = []; 
        centerSubregionContourSamples = 32;
        contourGenerationMethod = 'ellipseFitToPooledConeApertureImage';
        obj.generateVisualizationCache(xSupport, ySupport, centerSubregionContourSamples, contourGenerationMethod);
    end

    theContourData = obj.visualizationCache.rfCenterContourData{theRGCindex};

    % The cell position
    theCurrentRGCposition = obj.rgcRFpositionsDegs(theRGCindex,:);

    % Retrieve this cell's  center cone indices and weights
    connectivityVector = full(squeeze(obj.rgcRFcenterConePoolingMatrix(:, theRGCindex)));
    centerConeIndicesForCurrentRGC = find(connectivityVector > 0.0001);
    centerConeWeightsForCurrentRGC = connectivityVector(centerConeIndicesForCurrentRGC);

    % Retrieve this cell's surround cone indices and weights
    connectivityVector = full(squeeze(obj.rgcRFsurroundConePoolingMatrix(:, theRGCindex)));
    surroundConeIndicesForCurrentRGC = find(connectivityVector > 0.0001);
    surroundConeWeightsForCurrentRGC = connectivityVector(surroundConeIndicesForCurrentRGC);

    inputConeMosaic = obj.inputConeMosaic;

    spatialSupportRangeArcMin = tickSeparationArcMin*4;

    % Generate the axes
    if (isempty(theAxes))
        hFig = figure(); clf;
        ff = MSreadyPlot.figureFormat('1x3 RF poster');
        theAxes = MSreadyPlot.generateAxes(hFig,ff);
    end

    centerLineWeightingFunctions = mRGCMosaic.renderSubregionConePoolingPlot(theAxes{1,1}, ...
            inputConeMosaic, ...
            theCurrentRGCposition, ...
            centerConeIndicesForCurrentRGC, ...
            centerConeWeightsForCurrentRGC, ...
            'flatTopSaturationLevel', 1.0, ...
            'withFigureFormat', ff, ...
            'spatialSupportRangeArcMin', spatialSupportRangeArcMin, ...
            'tickSeparationArcMin', tickSeparationArcMin, ...
            'plotTitle', sprintf('RGCindex: %d', theRGCindex), ...
            'noXLabel', true, ...
            'noYLabel', true, ...
            'noYTicks', true, ...
            'xAxisTickAngleRotationDegs', 0);

    cla(theAxes{1,1}, 'reset');

    % Plot the RFs
    surroundLineWeightingFunctions = mRGCMosaic.renderSubregionConePoolingPlot(theAxes{1,1}, ...
            inputConeMosaic, ...
            theCurrentRGCposition, ...
            surroundConeIndicesForCurrentRGC, ...
            surroundConeWeightsForCurrentRGC, ...
            'flatTopSaturationLevel', surroundSaturationLevel, ...
            'withFigureFormat', ff, ...
            'spatialSupportRangeArcMin', spatialSupportRangeArcMin, ...
            'tickSeparationArcMin', tickSeparationArcMin, ...
            'plotTitle', '', ...
            'noXLabel', false, ...
            'noYLabel', false, ...
            'noYTicks', false, ...
            'xAxisTickAngleRotationDegs', 0);
    if (reverseXDir)
        set(theAxes{1,1}, 'XDir', 'reverse');
    end

    if (iscell(theContourData))
        S.Vertices = theContourData{1}.vertices;
        S.Faces = theContourData{1}.faces;
    else
        S.Vertices = theContourData.vertices;
        S.Faces = theContourData.faces;
    end

    S.FaceVertexCData = [0.5 0.5 0.5];

    S.FaceColor = [1 1 1];
    S.EdgeColor = [1 1 0];
    S.FaceAlpha = 0.0;
    S.LineWidth = 1.5;
    patch(S, 'Parent', theAxes{1,1});
    S.LineStyle = ':';
    S.EdgeColor = [0 0 0];
    S.LineWidth = 1.5;
    patch(S, 'Parent', theAxes{1,1});
    hold(theAxes{1,1}, 'off');


    % Add the line weighting functions
    sensitivityRange(2) =  max([max(centerLineWeightingFunctions.x.amplitude(:)) max(centerLineWeightingFunctions.y.amplitude(:))]);
    sensitivityRange(1) = -normalizedPeakSurroundSensitivity*sensitivityRange(2);

    centerLineWeightingFunctions.x.amplitude = centerLineWeightingFunctions.x.amplitude / max(sensitivityRange);
    centerLineWeightingFunctions.y.amplitude = centerLineWeightingFunctions.y.amplitude / max(sensitivityRange);
    surroundLineWeightingFunctions.x.amplitude = surroundLineWeightingFunctions.x.amplitude / max(sensitivityRange);
    surroundLineWeightingFunctions.y.amplitude = surroundLineWeightingFunctions.y.amplitude / max(sensitivityRange);
    sensitivityRange = sensitivityRange / max(sensitivityRange);

    if (~isempty(theAxes{1,2}))
        cla(theAxes{1,2}, 'reset');
        mRGCMosaic.renderSubregionConePoolingLineWeightingFunctions(theAxes{1,2}, ...
                centerLineWeightingFunctions.x, surroundLineWeightingFunctions.x, ...
                sensitivityRange, 'x', ...
                'withFigureFormat', ff, ...
                'spatialSupportRangeArcMin', spatialSupportRangeArcMin, ...
                'tickSeparationArcMin', tickSeparationArcMin, ...
                'plotTitle', '', ...
                'noYLabel', noYLabelLineWeightingFunctions, ...
                'noYTicks', true, ...
                'gridless', gridlessLineWeightingFunctions, ...
                'xAxisTickAngleRotationDegs', 0);
        if (reverseXDir)
            set(theAxes{1,2}, 'XDir', 'reverse');
        end
    end

    if (~isempty(theAxes{1,3}))
        cla(theAxes{1,3}, 'reset');
        mRGCMosaic.renderSubregionConePoolingLineWeightingFunctions(theAxes{1,3}, ...
                centerLineWeightingFunctions.y, surroundLineWeightingFunctions.y, ...
                sensitivityRange, 'y', ...
                'withFigureFormat', ff, ...
                'spatialSupportRangeArcMin', spatialSupportRangeArcMin, ...
                'tickSeparationArcMin', tickSeparationArcMin, ...
                'plotTitle', '', ...
                'noYLabel', noYLabelLineWeightingFunctions, ...
                'noYTicks', true, ...
                'gridless', gridlessLineWeightingFunctions, ...
                'xAxisTickAngleRotationDegs', 0);
    
        if (reverseXDir)
            set(theAxes{1,3}, 'XDir', 'reverse');
        end
    end

end