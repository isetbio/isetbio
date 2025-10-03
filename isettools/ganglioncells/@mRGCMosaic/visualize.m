function [hFig,ax] = visualize(obj, varargin)
    
    % Parse optional input
    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    
    p.addParameter('clearAxesBeforeDrawing', true, @islogical);
    p.addParameter('labelRetinalMeridians', false, @islogical);
    p.addParameter('domainVisualizationLimits', [], @(x)((isempty(x))||(numel(x)==4)));
    p.addParameter('domainVisualizationTicks', [], @(x)(isempty(x)||(isstruct(x)&&((isfield(x, 'x'))&&(isfield(x,'y'))))));
    p.addParameter('tickSeparationArcMinForRFconePoolingMap', 3, @isscalar);
    p.addParameter('visualizedRGCindices', 'all', @(x)(isnumeric(x)||(ischar(x)&&(strcmp(x,'all')))));
    p.addParameter('singleRGCconePoolingRFmaps', false, @islogical);
    p.addParameter('scaleBarDegs', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('doNotLabelScaleBar', false, @islogical);

    p.addParameter('centerSubregionContourSamples', 10, @isscalar);
    p.addParameter('spatialSupportSamples', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('forceRegenerateVisualizationCache', false, @islogical);
    p.addParameter('contourGenerationMethod', 'ellipseFitToPooledConeApertureImage', @(x)(ismember(x, obj.validRFsubregionContourGenerationMethods)));
    p.addParameter('minConeWeightVisualized', mRGCMosaic.sensitivityAtPointOfOverlap, @isscalar);
    p.addParameter('maxNumberOfConesOutsideContour', 1, @isscalar);
    p.addParameter('identifiedConeAperture', 'geometricArea', @(x)ismember(x, ...
        {'lightCollectingArea', 'geometricArea', 'coneSpacing', ...
        'lightCollectingAreaCharacteristicDiameter', 'lightCollectingArea2sigma', 'lightCollectingArea4sigma', 'lightCollectingArea5sigma', 'lightCollectingArea6sigma'}));
    p.addParameter('identifiedConeApertureThetaSamples', 8, @(x)(isempty(x) || isscalar(x)));
    p.addParameter('identifyInputCones', false, @islogical);
    p.addParameter('identifyPooledCones', false, @islogical);
    p.addParameter('identifiedInputConeIndices', [], @isnumeric);
    p.addParameter('identifiedInputConeIndicesContour', false, @islogical);
    p.addParameter('pooledConesLineWidth', 1.0, @(x)(isempty(x) || isscalar(x)));
    p.addParameter('inputConesAlpha', 0.7, @isscalar);
    p.addParameter('plotRFoutlines', true, @islogical);
    p.addParameter('plottedRFoutlineLineWidth', 1.0, @isscalar);
    p.addParameter('plottedRFoutlineFaceColor',  [], @(x)(isempty(x)||((isvector(x))&&(numel(x) == 3))));
    p.addParameter('plottedRFoutlineFaceAlpha', 1.0, @isscalar);
    p.addParameter('renderSourceLatticeInsteadOfConnectedRFcenters', false, @islogical);
    p.addParameter('withSuperimposedPSF', [], @(x)(isempty(x) || isstruct(x)));
    p.addParameter('activation', []);
    p.addParameter('activationRange', [],@(x)((isempty(x))||(numel(x)==2)));
    p.addParameter('activationColorMap', [], @(x)(isempty(x)||(size(x,2) == 3)));
    p.addParameter('horizontalActivationColorBar', false, @islogical);
    p.addParameter('verticalActivationColorBar', false, @islogical);
    p.addParameter('horizontalActivationColorBarInside', false, @islogical);
    p.addParameter('verticalActivationColorBarInside', false, @islogical);

    p.addParameter('fontSize', 16, @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('fontAngle', 'normal', @(x)(ismember(lower(x), {'normal', 'italic'})));
    p.addParameter('labelRGCsWithIndices', [], @(x)(isempty(x)||isnumeric(x)));
    p.addParameter('labeledRGCsColor', [], @(x)(isempty(x)||(size(x,2) == 3)));
    p.addParameter('labeledRGCsLineWidth', 1.5, @isscalar);
    p.addParameter('colorbarFontSize', 16, @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('colorBarTickLabelPostFix', '', @ischar);
    p.addParameter('colorbarTickLabelColor',  [], @(x)(isempty(x)||((isvector(x))&&(numel(x) == 3))));

    p.addParameter('backgroundColor', [], @(x)( (ischar(x)&&((strcmp(x,'none'))||(strcmp(x,'mean of color map'))) ) || isempty(x) || ((isvector(x))&&(numel(x) == 3))));
    p.addParameter('plotTitle', '', @(x)(isempty(x) || ischar(x) || islogical(x)));
    p.addParameter('plotTitleColor', [0 0 0], @isnumeric);
    p.addParameter('plotTitleFontSize', 16, @isscalar);
    
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);

    p.addParameter('superimposedRect', [], @(x)(isempty(x)||isstruct(x)));
    p.addParameter('superimposedRectLineWidth', 2, @isscalar);
    p.addParameter('superimposedRectColor', [], @(x)( isempty(x) || ((isvector(x))&&(numel(x) == 3))));
    p.addParameter('superimposedRectAlpha', 1, @(x)(isscalar(x)&&(x>=0)&&(x<=1)));

    p.addParameter('exportVisualizationPDF', false, @islogical);
    p.addParameter('exportVisualizationPDFdirectory', 'mosaicComponentVisualizations', @ischar);
    p.addParameter('visualizationPDFfileName', 'mRGCmosaic', @ischar);

    p.parse(varargin{:});

    hFig = p.Results.figureHandle;
    ax = p.Results.axesHandle;
    if (isempty(hFig)) && (~isempty(ax))
        try
            hFig = ancestor(ax, 'figure');
        catch
            fprintf('MRGCMosaic.visualize:: Unable to retrieve the parent figure handle from passed axes\n');
        end
    end

    visualizationPDFfileName = p.Results.visualizationPDFfileName;
    exportVisualizationPDF = p.Results.exportVisualizationPDF;
    exportVisualizationPDFdirectory = p.Results.exportVisualizationPDFdirectory;
    clearAxesBeforeDrawing = p.Results.clearAxesBeforeDrawing;
    labelRetinalMeridians = p.Results.labelRetinalMeridians;
    domainVisualizationLimits = p.Results.domainVisualizationLimits;
    domainVisualizationTicks = p.Results.domainVisualizationTicks;
    tickSeparationArcMinForRFconePoolingMap = p.Results.tickSeparationArcMinForRFconePoolingMap;

    if (~strcmp(p.Results.visualizedRGCindices, 'all'))
        visualizedRGCindices = unique(p.Results.visualizedRGCindices);
    else
        visualizedRGCindices = p.Results.visualizedRGCindices;
    end

    singleRGCconePoolingRFmaps = p.Results.singleRGCconePoolingRFmaps;
    scaleBarDegs = p.Results.scaleBarDegs;
    doNotLabelScaleBar = p.Results.doNotLabelScaleBar;

    spatialSupportSamples = p.Results.spatialSupportSamples;
    centerSubregionContourSamples = p.Results.centerSubregionContourSamples;
    contourGenerationMethod = p.Results.contourGenerationMethod;
    
    forceRegenerateVisualizationCache = p.Results.forceRegenerateVisualizationCache;
    minConeWeightVisualized = p.Results.minConeWeightVisualized;
    maxNumberOfConesOutsideContour = p.Results.maxNumberOfConesOutsideContour;
    identifiedConeAperture = p.Results.identifiedConeAperture;
    identifiedConeApertureThetaSamples = p.Results.identifiedConeApertureThetaSamples;
    identifiedInputConeIndices = p.Results.identifiedInputConeIndices;
    identifiedInputConeIndicesContour = p.Results.identifiedInputConeIndicesContour;
    identifyInputCones = p.Results.identifyInputCones;
    identifyPooledCones = p.Results.identifyPooledCones;
    pooledConesLineWidth = p.Results.pooledConesLineWidth;
    inputConesAlpha = p.Results.inputConesAlpha;
    plotRFoutlines = p.Results.plotRFoutlines;
    plottedRFoutlineLineWidth = p.Results.plottedRFoutlineLineWidth;
    plottedRFoutlineFaceColor = p.Results.plottedRFoutlineFaceColor;
    plottedRFoutlineFaceAlpha = p.Results.plottedRFoutlineFaceAlpha;

    renderSourceLatticeInsteadOfConnectedRFcenters = p.Results.renderSourceLatticeInsteadOfConnectedRFcenters;
    activation = p.Results.activation;
    activationRange = p.Results.activationRange;
    activationColorMap = p.Results.activationColorMap;
    colorBarTickLabelPostFix = p.Results.colorBarTickLabelPostFix;
    verticalColorBar = p.Results.verticalActivationColorBar;
    horizontalColorBar = p.Results.horizontalActivationColorBar;
    verticalColorBarInside = p.Results.verticalActivationColorBarInside;
    horizontalColorBarInside = p.Results.horizontalActivationColorBarInside;
    colorbarTickLabelColor = p.Results.colorbarTickLabelColor;
    backgroundColor = p.Results.backgroundColor;
    fontSize = p.Results.fontSize;
    fontAngle = p.Results.fontAngle;
    labelRGCsWithIndices = unique(p.Results.labelRGCsWithIndices);
    labeledRGCsColor = p.Results.labeledRGCsColor;
    labeledRGCsLineWidth = p.Results.labeledRGCsLineWidth;
    colorbarFontSize = p.Results.colorbarFontSize;
    plotTitle = p.Results.plotTitle;
    plotTitleColor = p.Results.plotTitleColor;
    plotTitleFontSize = p.Results.plotTitleFontSize;

    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;

    superimposedPSF = p.Results.withSuperimposedPSF;
    superimposedRect = p.Results.superimposedRect;
    superimposedRectLineWidth = p.Results.superimposedRectLineWidth;
    superimposedRectColor = p.Results.superimposedRectColor;
    superimposedRectAlpha = p.Results.superimposedRectAlpha;

    % validateSuperimposedRect
    if (~isempty(superimposedRect))
        assert(size(superimposedRect.center,2) == 2, 'superimposedRect must have 2 columns'); 
        assert(size(superimposedRect.center,1) == numel(superimposedRect.xRange), 'The rows of superimposedRect.center must equal the length of superimposedRect.xRange'); 
        assert(size(superimposedRect.center,1) == numel(superimposedRect.yRange), 'The rows of superimposedRect.center must equal the length of superimposedRect.xRange'); 
    end


    % Determine X,Y limits
    if (isempty(domainVisualizationLimits))
        if (isfield(obj.visualizationCache, 'surroundConesXrange'))
            xRange = obj.visualizationCache.surroundConesXrange;
        else
            xRange(1) = min(squeeze(obj.inputConeMosaic.coneRFpositionsDegs(:,1)));
            xRange(2) = max(squeeze(obj.inputConeMosaic.coneRFpositionsDegs(:,1)));
        end

        if (xRange(2) == xRange(1))
            xRange = xRange(1) + 0.02*[-1 1];
        end
        if (isfield(obj.visualizationCache, 'surroundConesYrange'))
            yRange = obj.visualizationCache.surroundConesYrange;
        else
            yRange(1) = min(squeeze(obj.inputConeMosaic.coneRFpositionsDegs(:,2)));
            yRange(2) = max(squeeze(obj.inputConeMosaic.coneRFpositionsDegs(:,2)));
        end

        if (yRange(2) == yRange(1))
            yRange = yRange(1) + 0.02*[-1 1];
        end
        xx = xRange(2)-xRange(1);
        yy = yRange(2)-yRange(1);
        XLims(1) = xRange(1)-xx*0.02;
        XLims(2) = xRange(2)+xx*0.02;
        YLims(1) = yRange(1)-yy*0.02;
        YLims(2) = yRange(2)+yy*0.02;
        domainVisualizationLimits(1) = XLims(1);
        domainVisualizationLimits(2) = XLims(2);
        domainVisualizationLimits(3) = YLims(1);
        domainVisualizationLimits(4) = YLims(2);
    else
        XLims(1) = domainVisualizationLimits(1);
        XLims(2) = domainVisualizationLimits(2);
        YLims(1) = domainVisualizationLimits(3);
        YLims(2) = domainVisualizationLimits(4);
    end
    
    if (isempty(domainVisualizationTicks))
        xo = (XLims(1)+XLims(2))/2;
        xx = XLims(2)-XLims(1);
        yo = (YLims(1)+YLims(2))/2;
        yy = YLims(2)-YLims(1);
        ticksX = [XLims(1) xo XLims(2)];
        ticksY = [YLims(1) yo YLims(2)];
        
        if (xx > 10)
            domainVisualizationTicks.x = round(ticksX);
        elseif (xx > 5)
            domainVisualizationTicks.x = round(ticksX*10)/10;
        elseif (xx > 1)
            domainVisualizationTicks.x = round(ticksX*100)/100;
        else
            domainVisualizationTicks.x = round(ticksX*1000)/1000;
        end
        if (yy > 10)
            domainVisualizationTicks.y = round(ticksY);
        elseif (yy > 5)
            domainVisualizationTicks.y = round(ticksY*10)/10;
        elseif (yy > 1)
            domainVisualizationTicks.y = round(ticksY*100)/100;
        else
            domainVisualizationTicks.y = round(ticksY*1000)/1000;
        end 
    end


    assert(domainVisualizationLimits(2)>domainVisualizationLimits(1), 'xRange: domainVisualizationLimits(2) must be > domainVisualizationLimits(1)');
    assert(domainVisualizationLimits(4)>domainVisualizationLimits(3), 'yRange: domainVisualizationLimits(4) must be > domainVisualizationLimits(2)');

    if (singleRGCconePoolingRFmaps)
        for iRGC = 1:numel(visualizedRGCindices)
            theRGCindex = visualizedRGCindices(iRGC);
            hFig = obj.renderConePoolingRFmap(theRGCindex, ...
                'figureHandle', hFig, ...
                'axesHandle', ax, ...
                'scaleBarDegs', scaleBarDegs, ...
                'doNotLabelScaleBar', doNotLabelScaleBar, ...
                'minConeWeightIncluded', minConeWeightVisualized, ...
                'tickSeparationArcMin', tickSeparationArcMinForRFconePoolingMap);

            if (exportVisualizationPDF)
                p = getpref('isetbio');
                pdfExportSubDir = fullfile(p.rgcResources.figurePDFsDir);
                theVisualizationPDFfilename = fullfile(exportVisualizationPDFdirectory, sprintf('RFmap%d.pdf', theRGCindex));
            
                % Generate the path if we need to
                RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                    pdfExportSubDir, theVisualizationPDFfilename, ...
                    'generateMissingSubDirs', true);
        
                thePDFfileName = fullfile(pdfExportSubDir, theVisualizationPDFfilename);
                NicePlot.exportFigToPDF(thePDFfileName, hFig, 300);
            end

        end
        return;
    end

    if (renderSourceLatticeInsteadOfConnectedRFcenters)
        contourGenerationMethod = 'ellipseFitBasedOnLocalSpacingOfSourceLattice';
        % Load data from the source lattice
        [~, obj.rgcRFpositionsDegsOfSourceLattice] = retinalattice.import.finalMRGCPositions(...
            obj.inputConeMosaic.sourceLatticeSizeDegs, ...
            obj.inputConeMosaic.eccentricityMicrons, ... 
            max(obj.inputConeMosaic.coneRFpositionsMicrons,[], 1)-min(obj.inputConeMosaic.coneRFpositionsMicrons,[], 1), ...
            obj.inputConeMosaic.whichEye, ...
            obj.inputConeMosaic.customMMsToDegsConversionFunction);
        visualizedRGCindices = 1:size(obj.rgcRFpositionsDegsOfSourceLattice,1);
        % Compute spacings
        obj.rgcRFspacingsDegsOfSourceLattice = RGCmodels.Watson.convert.positionsToSpacings(obj.rgcRFpositionsDegsOfSourceLattice);
        identifyPooledCones = false;
    end

    % Generate the visualization cache
    
    if (ischar(visualizedRGCindices))&&(strcmp(visualizedRGCindices,'all'))
        visualizedRGCindices = 1:obj.rgcsNum;
    end
    xSupport = []; ySupport = []; 
    obj.generateVisualizationCache(xSupport, ySupport, ...
        centerSubregionContourSamples, contourGenerationMethod, ...
        visualizedRGCindices, minConeWeightVisualized, ...
        'spatialSupportSamples', spatialSupportSamples, ...
        'maxNumberOfConesOutsideContour', maxNumberOfConesOutsideContour, ...
        'forceRegenerateVisualizationCache', forceRegenerateVisualizationCache);

    [hFig, ax] = renderMosaicOfRFcenters(obj, hFig, ax, clearAxesBeforeDrawing, ...
        labelRetinalMeridians, ...
        domainVisualizationTicks, domainVisualizationLimits, ...
        identifiedConeAperture, identifiedConeApertureThetaSamples, ...
        identifyInputCones, identifyPooledCones, identifiedInputConeIndices, identifiedInputConeIndicesContour, ...
        pooledConesLineWidth, inputConesAlpha, visualizedRGCindices, labelRGCsWithIndices, ...
        labeledRGCsColor, labeledRGCsLineWidth, ...
        scaleBarDegs, doNotLabelScaleBar, ...
        plotRFoutlines, plottedRFoutlineLineWidth, plottedRFoutlineFaceColor, plottedRFoutlineFaceAlpha, ...
        superimposedPSF, ...
        superimposedRect, superimposedRectLineWidth, superimposedRectColor, superimposedRectAlpha, ...
        activation, activationRange, activationColorMap, ...
        colorBarTickLabelPostFix, colorbarTickLabelColor, ...
        verticalColorBar, horizontalColorBar, colorbarFontSize, ...
        verticalColorBarInside, horizontalColorBarInside, ...
        noXLabel, noYLabel, ...
        backgroundColor, fontSize, fontAngle, ...
        plotTitle, plotTitleColor, plotTitleFontSize);


    if (exportVisualizationPDF)
        p = getpref('isetbio');
        pdfExportSubDir = fullfile(p.rgcResources.figurePDFsDir);
        theVisualizationPDFfilename = fullfile(exportVisualizationPDFdirectory, sprintf('%s.pdf', visualizationPDFfileName));
    
        % Generate the path if we need to
        RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
            pdfExportSubDir, theVisualizationPDFfilename, ...
            'generateMissingSubDirs', true);

        thePDFfileName = fullfile(pdfExportSubDir, theVisualizationPDFfilename);
        NicePlot.exportFigToPDF(thePDFfileName, hFig, 300);
    end

end


function [hFig, ax] = renderMosaicOfRFcenters(obj,hFig, ax, clearAxesBeforeDrawing, ...
        labelRetinalMeridians, ...
        domainVisualizationTicks, domainVisualizationLimits, ...
        identifiedConeAperture, identifiedConeApertureThetaSamples, ...
        identifyInputCones, identifyPooledCones, identifiedInputConeIndices, identifiedInputConeIndicesContour, ...
        pooledConesLineWidth, inputConesAlpha, visualizedRGCindices, labelRGCsWithIndices, ...
        labeledRGCsColor, labeledRGCsLineWidth, ...
        scaleBarDegs, doNotLabelScaleBar, ...
        plotRFoutlines, plottedRFoutlineLineWidth, plottedRFoutlineFaceColor, plottedRFoutlineFaceAlpha, ...
        superimposedPSF, ...
        superimposedRect, superimposedRectLineWidth, superimposedRectColor, superimposedRectAlpha, ...
        activation, activationRange, activationColorMap, ...
        colorBarTickLabelPostFix, colorbarTickLabelColor, ...
        verticalColorBar, horizontalColorBar, colorbarFontSize, ...
        verticalColorBarInside, horizontalColorBarInside, ...
        noXLabel, noYLabel, ...
        backgroundColor, fontSize, fontAngle, plotTitle, plotTitleColor,  plotTitleFontSize)


    if (isempty(ax))
        if (isempty(hFig))
            hFig = figure(); clf;
            set(hFig, 'Color', [1 1 1], 'Position', [10 10 1300 900]);
        end
        ax = subplot('Position', [0.05 0.06 0.91 0.91]);
    else
        if (clearAxesBeforeDrawing)
            cla(ax);
        end
    end

    if (~isempty(activation))
        % Visualize RF centers, color-coded with their activation level
        if (isempty(activationColorMap))
            linearRamp = linspace(0,1,1024);
            linearRamp = linearRamp(:);
            cMap = [linearRamp*0 linearRamp linearRamp*0];
        else
            cMap = activationColorMap;
        end
        
        if (ischar(backgroundColor) && strcmp(backgroundColor, 'mean of color map'))
            midRow = round(size(cMap,1)/2);
            backgroundColor = squeeze(cMap(midRow,:));
        elseif (isempty(backgroundColor))
            backgroundColor = squeeze(cMap(1,:));
        end
        
        if (isempty(activationRange))
            activationRange = [min(activation(:)) max(activation(:))];
        end

        activation = (activation - activationRange(1))/(activationRange(2)-activationRange(1));
        activation(activation<0) = 0;
        activation(activation>1) = 1;

        currentFacesNum = 0;
        for iRGC = 1:obj.rgcsNum
            newVerticesNum = obj.visualizationCache.rfCenterPatchData.verticesNumForRGC(iRGC);
            idx = currentFacesNum + (1:newVerticesNum);
            obj.visualizationCache.rfCenterPatchData.faceVertexCData(idx,:) = activation(iRGC);
            currentFacesNum = currentFacesNum + newVerticesNum;
        end
    else
        if (isempty(backgroundColor))
            backgroundColor = [1 1 1];
        end

        % All mRGC centers in gray
        cMap = [0 0 0; 0.5 0.5 0.5; 0 0 0];
        S.FaceVertexCData = obj.visualizationCache.rfCenterPatchData.faceVertexCData*0+0.5;
    end

    if (plotRFoutlines) || (~isempty(activation))
        % Plot the RFs
        S.Vertices = obj.visualizationCache.rfCenterPatchData.vertices;
        S.Faces = obj.visualizationCache.rfCenterPatchData.faces;
        S.FaceVertexCData = obj.visualizationCache.rfCenterPatchData.faceVertexCData;
    
        S.FaceColor = 'flat';
        if (~isempty(activation))
            S.EdgeColor = [0 0 0]';
        else
            if (isempty(plottedRFoutlineFaceColor))
                S.FaceColor = [0.95 0.95 0.95]*0.6;
                S.EdgeColor = [0 0 0];
            else
                S.FaceColor = plottedRFoutlineFaceColor;
                S.EdgeColor = S.FaceColor*0.25;
            end
        end

        if (isempty(activation))
            S.FaceAlpha = plottedRFoutlineFaceAlpha;
        else
            S.FaceAlpha = 1.0;
        end

        S.EdgeAlpha = 1.0;
        S.LineWidth = plottedRFoutlineLineWidth;
        patch(S, 'Parent', ax);
    end

    if (identifyPooledCones)
        hold(ax, 'on')
        lConeInputLineColor = [1 0.1 0.5];
        mConeInputLineColor = [0.1 1 0.5];

        % Put a single dot in all mRGC RF centers with a single input
        if (~isempty(obj.visualizationCache.rfCenterSingleConeInputDotPositions))
            plot(ax, obj.visualizationCache.rfCenterSingleConeInputDotPositions(:,1), ...
                     obj.visualizationCache.rfCenterSingleConeInputDotPositions(:,2), 'k.');
        end
     
        % Render line segments from centroid to pulled cones
        renderPooledConesLineSegments(obj, ax, lConeInputLineColor, mConeInputLineColor, pooledConesLineWidth);
    end

    % Identify input cones
    if (identifyInputCones)
        hold(ax, 'on')
        if (isempty(identifiedInputConeIndices))
            if (identifiedConeApertureThetaSamples > 0)
                obj.inputConeMosaic.visualize(...
                    'figureHandle', hFig, 'axesHandle', ax, ...
                    'clearAxesBeforeDrawing', false, ...
                    'visualizedConeAperture', identifiedConeAperture, ...
                    'conesAlpha', inputConesAlpha, ...
                    'conesEdgeAlpha', inputConesAlpha, ...
                    'visualizedConeApertureThetaSamples', identifiedConeApertureThetaSamples, ...
                    'withSuperimposedPSF', superimposedPSF, ...
                    'labelRetinalMeridians', labelRetinalMeridians, ...
                    'domainVisualizationTicks', domainVisualizationTicks, ...
                    'domainVisualizationLimits', domainVisualizationLimits, ...
                    'backgroundColor', backgroundColor);
            else
                plot(ax, obj.inputConeMosaic.coneRFpositionsDegs(obj.inputConeMosaic.lConeIndices,1), ...
                         obj.inputConeMosaic.coneRFpositionsDegs(obj.inputConeMosaic.lConeIndices,2), ...
                         'r.', 'MarkerSize', 16, 'LineWidth', 1.0);
                plot(ax, obj.inputConeMosaic.coneRFpositionsDegs(obj.inputConeMosaic.mConeIndices,1), ...
                         obj.inputConeMosaic.coneRFpositionsDegs(obj.inputConeMosaic.mConeIndices,2), ...
                         'g.', 'MarkerSize', 16,'LineWidth', 1.0);
                plot(ax, obj.inputConeMosaic.coneRFpositionsDegs(obj.inputConeMosaic.sConeIndices,1), ...
                         obj.inputConeMosaic.coneRFpositionsDegs(obj.inputConeMosaic.sConeIndices,2), ...
                         'b.', 'MarkerSize', 16, 'LineWidth', 1.0);
                set(ax, 'XLim', [domainVisualizationLimits(1) domainVisualizationLimits(2)], ...
                        'YLim', [domainVisualizationLimits(3) domainVisualizationLimits(4)]);
                set(ax, 'XTick', domainVisualizationTicks.x, 'YTick', domainVisualizationTicks.y);
            end
        else
            obj.inputConeMosaic.visualize(...
                'figureHandle', hFig, 'axesHandle', ax, ...
                'clearAxesBeforeDrawing', false, ...
                'visualizedConeAperture', identifiedConeAperture, ...
                'labelConesWithIndices', identifiedInputConeIndices, ...
                'conesAlpha', inputConesAlpha, ...
                'conesEdgeAlpha', inputConesAlpha, ...
                'visualizedConeApertureThetaSamples', identifiedConeApertureThetaSamples, ...
                'withSuperimposedPSF', superimposedPSF, ...
                'labelRetinalMeridians', labelRetinalMeridians, ...
                'domainVisualizationTicks', domainVisualizationTicks, ...
                'domainVisualizationLimits', domainVisualizationLimits, ...
                'backgroundColor', backgroundColor);

            if (identifiedInputConeIndicesContour)
                hold(ax, 'on')
                if (numel(visualizedRGCindices) == 1)
                    % Add dashed line showing the surround extent
                    theRGCpos = obj.rgcRFpositionsDegs(visualizedRGCindices,:);
                    maxRadius = max(sqrt(sum(bsxfun(@minus,obj.inputConeMosaic.coneRFpositionsDegs(identifiedInputConeIndices,:), theRGCpos).^2,2)));
                    xx = theRGCpos(1) + maxRadius*cosd(0:10:350);
                    yy = theRGCpos(2) + maxRadius*sind(0:10:350);
                    plot(ax, xx,yy, 'k--', 'LineWidth', 2.0);
                end
            end
        end
    else
        hold(ax, 'on')
        obj.inputConeMosaic.visualize(...
            'figureHandle', hFig, 'axesHandle', ax, ...
            'clearAxesBeforeDrawing', false, ...
            'labelCones', false, ...
            'visualizedConeApertureThetaSamples', identifiedConeApertureThetaSamples, ...
            'withSuperimposedPSF', superimposedPSF, ...
            'labelRetinalMeridians', labelRetinalMeridians, ...
            'domainVisualizationTicks', domainVisualizationTicks, ...
            'domainVisualizationLimits', domainVisualizationLimits, ...
            'backgroundColor', backgroundColor);
    end


    if (~isempty(labelRGCsWithIndices))
        if (plotRFoutlines) || (~isempty(activation))
            hold(ax, 'on')
            for iRGC = 1:numel(labelRGCsWithIndices)
                theLabeledRGCindex = find(visualizedRGCindices == labelRGCsWithIndices(iRGC));
                if (~isempty(theLabeledRGCindex))
                    if (iscell(obj.visualizationCache.rfCenterContourData{theLabeledRGCindex}))
                        S = obj.visualizationCache.rfCenterContourData{theLabeledRGCindex}{1};
                    else
                        S = obj.visualizationCache.rfCenterContourData{theLabeledRGCindex};
                    end
    
                    S.FaceVertexCData = 0.5;
                    S.FaceColor = 'flat';
                    
    
                    S.FaceAlpha = 0.0;
                    S.LineWidth = 2*labeledRGCsLineWidth;
                    
                    S.EdgeColor = [0 0 0];
                    S.LineStyle = '-';
                    patch(S, 'Parent', ax);

                    if (isempty(labeledRGCsColor))
                        S.EdgeColor = [1 1 0];
                    else
                        S.EdgeColor = labeledRGCsColor;
                    end

                    S.LineWidth = labeledRGCsLineWidth;
                    S.LineStyle = '--';
                    patch(S, 'Parent', ax);
                end
            end
        end
    end


    if (~isempty(superimposedRect))
        hold(ax, 'on');
        if (isempty(superimposedRectColor))
            superimposedRectColor = [1 0 0];
        end

        for iRect = 1:size(superimposedRect.center,1)
            x1 = superimposedRect.center(iRect,1) - 0.5*superimposedRect.xRange(iRect);
            x2 = superimposedRect.center(iRect,1) + 0.5*superimposedRect.xRange(iRect);
            y1 = superimposedRect.center(iRect,2) - 0.5*superimposedRect.yRange(iRect);
            y2 = superimposedRect.center(iRect,2) + 0.5*superimposedRect.yRange(iRect);
            xx = [x1 x1 x2 x2 x1];
            yy = [y1 y2 y2 y1 y1];
            patchVertices = [xx(:) yy(:)];
            patchFaces = 1:size(patchVertices,1);
            patch('Faces', patchFaces, 'Vertices', patchVertices, ...
                'FaceColor', superimposedRectColor,  'FaceAlpha', superimposedRectAlpha, 'Parent', ax);

            plot(ax, xx,yy, '-', 'LineWidth', superimposedRectLineWidth, 'Color', ones(size(superimposedRectColor)) - superimposedRectColor);

            plot(ax, xx,yy, '--', 'LineWidth', superimposedRectLineWidth, 'Color', superimposedRectColor);
        end % iRect
    end

    % Add a scale bar for comparison with physiology
    if (~isempty(scaleBarDegs)) && (scaleBarDegs > 0) 
        hold(ax, 'on');
        xOffset = domainVisualizationLimits(1)+0.05*(domainVisualizationLimits(2)-domainVisualizationLimits(1));
        yOffset = domainVisualizationLimits(3)+0.05*(domainVisualizationLimits(4)-domainVisualizationLimits(3));
        yOffset2 = yOffset + 0.05*(domainVisualizationLimits(4)-domainVisualizationLimits(3));
        plot(ax, xOffset+[0 scaleBarDegs], yOffset2*[1 1], 'k-', 'LineWidth', 2);

        if (~doNotLabelScaleBar)
            if (scaleBarDegs>=1.0)
                text(ax, xOffset+scaleBarDegs, yOffset, sprintf(' %2.1f degs', scaleBarDegs), ...
                    'FontSize', 20, 'Color', [0 0 0], 'BackgroundColor', [0.85 0.85 0.85]);
            elseif (scaleBarDegs>=0.1)
                text(ax, xOffset+scaleBarDegs, yOffset, sprintf(' %2.2f degs', scaleBarDegs), ...
                    'FontSize', 20, 'Color', [0 0 0], 'BackgroundColor', [0.85 0.85 0.85]);
            else
                text(ax, xOffset+scaleBarDegs, yOffset, sprintf(' %2.3f degs', scaleBarDegs), ...
                    'FontSize', 20, 'Color', [0 0 0], 'BackgroundColor', [0.85 0.85 0.85]);
            end
        end
    end

    % Finalize plot
    hold(ax, 'off');
    set(ax, 'FontSize', fontSize, 'FontAngle', fontAngle);

    minTickIncrement = min([min(abs(diff(domainVisualizationTicks.x))) min(abs(diff(domainVisualizationTicks.y)))]);
    set(ax, 'XTick', domainVisualizationTicks.x, 'YTick', domainVisualizationTicks.y);
    if (minTickIncrement > 1-100*eps)
       set(ax, 'XTickLabel', sprintf('%1.0f\n', sign(domainVisualizationTicks.x).*round(abs(domainVisualizationTicks.x))), ...
               'YTickLabel', sprintf('%1.0f\n', sign(domainVisualizationTicks.y).*round(abs(domainVisualizationTicks.y))));
    elseif (minTickIncrement >= 0.1-100*eps)
       set(ax, 'XTickLabel', sprintf('%1.1f\n', domainVisualizationTicks.x), ...
               'YTickLabel', sprintf('%1.1f\n', domainVisualizationTicks.y));
    elseif (minTickIncrement > 0.01-100*eps)
       set(ax, 'XTickLabel', sprintf('%1.2f\n', domainVisualizationTicks.x), ...
               'YTickLabel', sprintf('%1.2f\n', domainVisualizationTicks.y));
    else
        set(ax, 'XTickLabel', sprintf('%1.3f\n', domainVisualizationTicks.x), ...
                'YTickLabel', sprintf('%1.3f\n', domainVisualizationTicks.y));
    end

    set(ax, 'XLim', domainVisualizationLimits(1:2), 'YLim', domainVisualizationLimits(3:4));
    if (~identifyInputCones)
        colormap(ax, cMap);
    end

    if (isempty(backgroundColor))
        set(ax, 'CLim', [0 1], 'Color', 'none');
    else
        set(ax, 'CLim', [0 1], 'Color', backgroundColor);
    end

    if (isempty(colorbarTickLabelColor))
        colorbarTickLabelColor = [1 0.5 0];
    end
    
    box(ax, 'on')

    % Colorbar and colorbar ticks
    if (~isempty(activation))
        if (verticalColorBar) || (horizontalColorBar) || (verticalColorBarInside) || (horizontalColorBarInside)
            colorBarTicks = [0.00 0.25 0.5 0.75 1.0];
            colorBarTickLabels = cell(1, numel(colorBarTicks));
            colorBarTickLevels = activationRange(1) + (activationRange(2)-activationRange(1)) * colorBarTicks;
            
            for k = 1:numel(colorBarTicks)
                if (max(abs(colorBarTickLevels)) >= 10)
                    colorBarTickLabels{k} = sprintf('%2.0f %s', colorBarTickLevels(k), colorBarTickLabelPostFix);
                elseif (max(abs(colorBarTickLevels)) >= 1)
                    colorBarTickLabels{k} = sprintf('%2.1f %s', colorBarTickLevels(k), colorBarTickLabelPostFix);
                elseif (max(abs(colorBarTickLevels)) > 0.1)
                    colorBarTickLabels{k} = sprintf('%2.2f %s', colorBarTickLevels(k), colorBarTickLabelPostFix);
                else
                    colorBarTickLabels{k} = sprintf('%2.3f %s', colorBarTickLevels(k), colorBarTickLabelPostFix);
                end
            end
            
            if (isempty(colorbarFontSize))
                colorbarFontSize = fontSize/2;
            end
    
  
            if (verticalColorBar)
                colorbar(ax, 'eastOutside', 'Ticks', colorBarTicks, 'TickLabels', colorBarTickLabels, ...
                    'Color', colorbarTickLabelColor);
            elseif (verticalColorBarInside)
                colorbar(ax, 'east', 'Ticks', colorBarTicks, 'TickLabels', colorBarTickLabels, ...
                    'Color', colorbarTickLabelColor,  'FontWeight', 'Bold', 'FontSize', colorbarFontSize, 'FontName', 'Spot mono');
            elseif (horizontalColorBar)
                colorbar(ax,'northOutside', 'Ticks', colorBarTicks, 'TickLabels', colorBarTickLabels, ...
                    'Color', colorbarTickLabelColor);
            elseif (horizontalColorBarInside)
                colorbar(ax,'north', 'Ticks', colorBarTicks, 'TickLabels', colorBarTickLabels, ...
                    'Color', colorbarTickLabelColor,  'FontWeight', 'Bold', 'FontSize', colorbarFontSize, 'FontName', 'Spot mono');
            end
        else
            colorbar(ax, 'off');
        end
    else
        colorbar(ax, 'off');
    end

    if (noXLabel)
        xlabel(ax, '');
    end
    if (noYLabel)
        ylabel(ax, '');
    end


    title(ax, plotTitle, 'Color', plotTitleColor, 'FontSize', plotTitleFontSize);
    drawnow;
end


function renderPooledConesLineSegments(obj,ax, lConeInputLineColor, mConeInputLineColor, lineSegmentWidth)
   
    if (~isempty(obj.visualizationCache.rfCenterConeConnectionLineSegments))
        % Plot the connections from the RF center to the input L-cones
        idx = find(obj.visualizationCache.rfCenterConeConnectionLineSegments.coneTypes == cMosaic.LCONE_ID);
    
        if (isempty(lineSegmentWidth))
            for i = 1:numel(idx)
                plot(ax, ...
                obj.visualizationCache.rfCenterConeConnectionLineSegments.Xpos(:, idx(i)), ...
                obj.visualizationCache.rfCenterConeConnectionLineSegments.Ypos(:,idx(i)), ...
                'Color', [0 0 0],...
                'LineWidth', obj.visualizationCache.rfCenterConeConnectionLineSegments.lineSegmentWidths(idx(i))*2); 

                plot(ax, ...
                obj.visualizationCache.rfCenterConeConnectionLineSegments.Xpos(:, idx(i)), ...
                obj.visualizationCache.rfCenterConeConnectionLineSegments.Ypos(:,idx(i)), ...
                'Color', lConeInputLineColor,...
                'LineWidth', obj.visualizationCache.rfCenterConeConnectionLineSegments.lineSegmentWidths(idx(i)));
            end
        else
            plot(ax, ...
                obj.visualizationCache.rfCenterConeConnectionLineSegments.Xpos(:, idx), ...
                obj.visualizationCache.rfCenterConeConnectionLineSegments.Ypos(:,idx), ...
                'Color', [0 0 0],...
                'LineWidth', lineSegmentWidth*2);  

            plot(ax, ...
                obj.visualizationCache.rfCenterConeConnectionLineSegments.Xpos(:, idx), ...
                obj.visualizationCache.rfCenterConeConnectionLineSegments.Ypos(:,idx), ...
                'Color', lConeInputLineColor,...
                'LineWidth', lineSegmentWidth); 
        end
    
        idx = find(obj.visualizationCache.rfCenterConeConnectionLineSegments.coneTypes == cMosaic.MCONE_ID);

        if (isempty(lineSegmentWidth))
            for i = 1:numel(idx)
                plot(ax, ...
                    obj.visualizationCache.rfCenterConeConnectionLineSegments.Xpos(:, idx(i)), ...
                    obj.visualizationCache.rfCenterConeConnectionLineSegments.Ypos(:,idx(i)), ...
                    'Color', [0 0 0],...
                    'LineWidth', obj.visualizationCache.rfCenterConeConnectionLineSegments.lineSegmentWidths(idx(i))*2);
                plot(ax, ...
                    obj.visualizationCache.rfCenterConeConnectionLineSegments.Xpos(:, idx(i)), ...
                    obj.visualizationCache.rfCenterConeConnectionLineSegments.Ypos(:,idx(i)), ...
                    'Color', mConeInputLineColor,...
                    'LineWidth', obj.visualizationCache.rfCenterConeConnectionLineSegments.lineSegmentWidths(idx(i)));
            end
        else
            plot(ax, ...
                obj.visualizationCache.rfCenterConeConnectionLineSegments.Xpos(:, idx), ...
                obj.visualizationCache.rfCenterConeConnectionLineSegments.Ypos(:,idx), ...
                'Color', [0 0 0],...
                'LineWidth', lineSegmentWidth*2);
            plot(ax, ...
                obj.visualizationCache.rfCenterConeConnectionLineSegments.Xpos(:, idx), ...
                obj.visualizationCache.rfCenterConeConnectionLineSegments.Ypos(:,idx), ...
                'Color', mConeInputLineColor,...
                'LineWidth', lineSegmentWidth);
        end
    end

end