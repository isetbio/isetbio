function [centerLineWeightingFunctions, surroundLineWeightingFunctions] = singleRGCconePoolingMap(figNo, theMRGCMosaic, theRGCindex, pdfFileName, varargin)
    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('figureFormat', [], @(x)(isempty(x)||isstruct(x)));
    p.addParameter('noGrid', false, @islogical);
    p.addParameter('minCenterConeWeight', 0.1, @isscalar);
    p.addParameter('minSurroundConeWeight', 0.5/100, @isscalar);
    p.addParameter('minSurroundConeWeightRelativity', 'center', @(x)(ismember(x, {'center', })));
    p.addParameter('domainVisualizationLimits', [], @(x)(isempty(x)||(numel(x)==4)));
    p.addParameter('domainVisualizationTicks', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('fixedSpatialSupportTickSeparationArcMin', @(x)(isempty(x)||isscalar(x)));
    p.addParameter('fixedScaleBarDegs', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('doNotLabelScaleBar', false, @islogical);
    p.addParameter('contourGenerationMethod', 'ellipseFitToPooledConeApertureImage', ...
        @(x)(ismember(x, mRGCMosaic.validRFsubregionContourGenerationMethods)));
    p.addParameter('maxNumberOfConesOutsideContour', 1, @isscalar);
    p.addParameter('plotTitle', '', @(x)(isempty(x)||ischar(x)));
    p.addParameter('noTitle', false, @islogical);
    p.addParameter('renderLineWeightingFunctionPlots', true, @islogical);
    p.parse(varargin{:});

    figureHandle = p.Results.figureHandle;
    axesHandle = p.Results.axesHandle;
    figureFormat = p.Results.figureFormat;
    noGrid = p.Results.noGrid;
    domainVisualizationLimits = p.Results.domainVisualizationLimits;
    domainVisualizationTicks = p.Results.domainVisualizationTicks;
    fixedSpatialSupportTickSeparationArcMin = p.Results.fixedSpatialSupportTickSeparationArcMin;
    fixedScaleBarDegs = p.Results.fixedScaleBarDegs;
    doNotLabelScaleBar = p.Results.doNotLabelScaleBar;

    maxNumberOfConesOutsideContour = p.Results.maxNumberOfConesOutsideContour;
    contourGenerationMethod = p.Results.contourGenerationMethod;
    renderLineWeightingFunctionPlots = p.Results.renderLineWeightingFunctionPlots;
    plotTitle = p.Results.plotTitle;
    noTitle = p.Results.noTitle;

    theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    pdfExportSubDir = 'validation';

    if (~isempty(figureHandle)) && (~isempty(axesHandle)) && (~isempty(figureFormat))
        hFig = figureHandle;
        ax = axesHandle;
        ff = figureFormat;
    else
        ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure'); 
        hFig = figure(figNo); clf;
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
        ax = theAxes{1,1};
    end


    % Criterion used by Field (2010) for center cones is 0.1:
    % (from supplemental methods section titled: "Cones providing input to receptive field center and surround"
    % "Cones in the RF center were defined as those with a weight greater than or equal to 10% of the weight of cone 
    % with the largest weight (which is by definition in the center)."
    % minCenterConeWeight = exp(-0.5);
    minCenterConeWeight = p.Results.minCenterConeWeight;

    % Criterion used by Field (2010) for surround cones is 0.5% of the center cone weight:
    % (from supplemental methods section titled: "Cones providing input to receptive field center and surround"
    % "Surround cones were defined as those cones with a light response polarity opposite to the center, 
    % weights >0.5% that of the peak cone but of opposite sign, and within 8 SDs of a circular Gaussian fit to the RF center
    % Note however, that their measured center weights are composite, center-surround as you cant
    % separate the surround response from the center response electrophysiologically, so we should do the same here)
    minSurroundConeWeightRelativity = p.Results.minSurroundConeWeightRelativity;
    minSurroundConeWeight = p.Results.minSurroundConeWeight


    % The spatial support
    spatialSupportCenterDegs = theMRGCMosaic.rgcRFpositionsDegs(theRGCindex,:);
    ecc = sqrt(sum(spatialSupportCenterDegs.^2,2));

    if (ecc > 20)
        scaleBarDegs = 0.1;
        spatialSupportTickSeparationArcMin = 30.0;
    elseif (ecc > 10)
        scaleBarDegs = 0.08;
        spatialSupportTickSeparationArcMin = 16.0;
    elseif (ecc > 5)
        scaleBarDegs = 0.07;
        spatialSupportTickSeparationArcMin = 12.0;
    elseif (ecc > 2)
        scaleBarDegs = 0.06;
        spatialSupportTickSeparationArcMin = 6.0;
    else
        scaleBarDegs = 0.05;
        spatialSupportTickSeparationArcMin = 6.0;
    end

    if (~isempty(fixedSpatialSupportTickSeparationArcMin))
        spatialSupportTickSeparationArcMin = fixedSpatialSupportTickSeparationArcMin;
    end

    if (~isempty(fixedScaleBarDegs))
        scaleBarDegs = fixedScaleBarDegs;
    end

    % surround alone (together with center ellipse)
    % visualizedConeWeights = 'center-surround';  % center - surround
    visualizedConeWeights = 'surround-alone';     % only the surround
    if (isempty(plotTitle))
        plotTitle = sprintf('visualized cone weights: ''%s''\n Wc > %1.3f, Ws > %1.3f (x peak %s cone)',...
            visualizedConeWeights, minCenterConeWeight, minSurroundConeWeight, minSurroundConeWeightRelativity);
    end

    if (noTitle)
        plotTitle = '';
    end

    theCenterConeIndices = theMRGCMosaic.singleCellConnectivityStats(...
        theRGCindex, 'center', ...
        'minConeWeightIncluded', 0.0001, ...
        'inputConeIndicesOnly', true, ...
        'warnIfCenterConeInputNumerosityDiffersFromExclusiveOne', false);
    theCenterConeWeights = full(squeeze(theMRGCMosaic.rgcRFcenterConeConnectivityMatrix(theCenterConeIndices, theRGCindex)));

    theSurroundConeIndices = theMRGCMosaic.singleCellConnectivityStats(...
        theRGCindex, 'surround', ...
        'minConeWeightIncluded', 0.0001, ...
        'inputConeIndicesOnly', true);
    theSurroundConeWeights = full(squeeze(theMRGCMosaic.rgcRFsurroundConeConnectivityMatrix(theSurroundConeIndices, theRGCindex)));

    thePooledConeIndicesAndWeightsStruct = struct();
    thePooledConeIndicesAndWeightsStruct.centerConeIndices = theCenterConeIndices;
    thePooledConeIndicesAndWeightsStruct.centerConeWeights = theCenterConeWeights;
    thePooledConeIndicesAndWeightsStruct.surroundConeIndices = theSurroundConeIndices;
    thePooledConeIndicesAndWeightsStruct.surroundConeWeights = theSurroundConeWeights;
    thePooledConeIndicesAndWeightsStruct.unconnectableSurroundConeIndices = [];
    thePooledConeIndicesAndWeightsStruct.unonnectableSurroundConeWeights = [];

    figPos = [];
    [centerLineWeightingFunctions, surroundLineWeightingFunctions] = ...
         RGCMosaicConstructor.visualize.conePoolingWeightsMap(...
            pdfExportSubDir, figNo, figPos, ...
            theMRGCMosaic, theRGCindex, ...
            spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
            thePooledConeIndicesAndWeightsStruct , ...
            minCenterConeWeight, ...
            minSurroundConeWeight, ...
            minSurroundConeWeightRelativity, ...
            visualizedConeWeights, ...
            'domainVisualizationLimits', domainVisualizationLimits, ...
            'domainVisualizationTicks', domainVisualizationTicks, ...
            'contourGenerationMethod', contourGenerationMethod, ...
            'maxNumberOfConesOutsideContour', maxNumberOfConesOutsideContour, ...
            'scaleBarDegs', scaleBarDegs, ...
            'doNotLabelScaleBar', doNotLabelScaleBar, ...
            'plotTitle', plotTitle, ...
            'axesToRenderIn', ax);

    ff.box = 'on';
    ff.tickDir = 'in';

    if (noGrid)
        ff.grid = 'off';
    end

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);

    if (~isempty(pdfFileName))
        fprintf('Empty pdf filename. Not exporting figure\n');
        thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, pdfFileName);
        NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
    end

    % Render line weighting function plots if so desired
    if (renderLineWeightingFunctionPlots)
        % Now the horizontal line weighting function
        whichMeridian = 'horizontal';
        theLineWeightingFunctionPDFfilename = strrep(...
            thePDFfileName, 'conePoolingWeightsMap', 'horizontalConePoolingLineWeightingFunction');

        hFig = figure(figNo+10); clf;
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
        ax = theAxes{1,1};
        if (~isempty(domainVisualizationLimits) && ~isempty(domainVisualizationTicks))
            RGCMosaicConstructor.visualize.centerAndSurroundConePoolingLineWeightingFunctions(...
                pdfExportSubDir, figNo, figPos, ...
                spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
                centerLineWeightingFunctions, surroundLineWeightingFunctions, whichMeridian, ...
                'axesToRenderIn', ax, ...
                'domainVisualizationLimits', domainVisualizationLimits(1:2), ...
                'domainVisualizationTicks', domainVisualizationTicks.x);
        else
            RGCMosaicConstructor.visualize.centerAndSurroundConePoolingLineWeightingFunctions(...
                pdfExportSubDir, figNo, figPos, ...
                spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
                centerLineWeightingFunctions, surroundLineWeightingFunctions, whichMeridian, ...
                'axesToRenderIn', ax);
        end

        ff.box = 'on';
        ff.tickDir = 'in';

        if (noGrid)
            ff.grid = 'off';
        end


        % Finalize figure using the Publication-Ready format
        PublicationReadyPlotLib.applyFormat(ax,ff);

        if (~isempty(pdfFileName))
            NicePlot.exportFigToPDF(theLineWeightingFunctionPDFfilename,hFig,  300);
        end

        % And the vertical line weighting function
        whichMeridian = 'vertical';
        theLineWeightingFunctionPDFfilename = strrep(...
            thePDFfileName, 'conePoolingWeightsMap', 'verticalConePoolingLineWeightingFunction');


        hFig = figure(figNo+20); clf;
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
        ax = theAxes{1,1};
        if (~isempty(domainVisualizationLimits) && ~isempty(domainVisualizationTicks))
            RGCMosaicConstructor.visualize.centerAndSurroundConePoolingLineWeightingFunctions(...
                pdfExportSubDir, figNo, figPos, ...
                spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
                centerLineWeightingFunctions, surroundLineWeightingFunctions, whichMeridian, ...
                'axesToRenderIn', ax, ...
                'noGrid', true, ...
                'domainVisualizationLimits', domainVisualizationLimits(3:4), ...
                'domainVisualizationTicks', domainVisualizationTicks.y);
        else
            RGCMosaicConstructor.visualize.centerAndSurroundConePoolingLineWeightingFunctions(...
                pdfExportSubDir, figNo, figPos, ...
                spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
                centerLineWeightingFunctions, surroundLineWeightingFunctions, whichMeridian, ...
                'axesToRenderIn', ax, ...
                'noGrid', true);
        end

        if (noGrid)
            ff.grid = 'off';
        end
        
        % Finalize figure using the Publication-Ready format
        PublicationReadyPlotLib.applyFormat(ax,ff);

        if (~isempty(pdfFileName))
            NicePlot.exportFigToPDF(theLineWeightingFunctionPDFfilename,hFig,  300);
        end
    end

end
