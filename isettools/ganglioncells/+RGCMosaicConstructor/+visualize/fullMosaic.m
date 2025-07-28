function [hFig, theAxes] = fullMosaic(theCenterConnectedMRGCMosaicFullFileName, theCenterConnectedMRGCMosaicFileName, ...
        minConeWeightVisualized, pStruct, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('contourGenerationMethod', 'ellipseFitToPooledConePositions', @(x)(ismember(x, mRGCMosaic.validRFsubregionContourGenerationMethods)));
    p.addParameter('visualizeRGCsAroundPositionAndWithingPositionDegs', [], @(x)(isempty(x)||(isnumeric(x)&&(numel(x)==4))));
    p.addParameter('identifyInputCones', false, @islogical);
    p.addParameter('identifyPooledCones', false, @islogical);
    p.addParameter('domainVisualizationLimits', [], @(x)(isempty(x)||(isnumeric(x)&&(numel(x)==4))));
    p.addParameter('domainVisualizationTicks', [], @(x)(isempty(x)||isstruct(x)));
    p.addParameter('visualizedSamplingPositionsGrid', [], @(x)(isempty(x)||isnumeric(x)));
    p.addParameter('visualizedSamplingPositionsGridColor', [1 1 0], @(x)(isnumeric(x)&&(size(x,2)==3)));
    p.addParameter('withFigureFormat', [], @isstruct);
    p.parse(varargin{:});
    visualizeRGCsAroundPositionAndWithingPositionDegs = p.Results.visualizeRGCsAroundPositionAndWithingPositionDegs;
    identifyInputCones = p.Results.identifyInputCones;
    identifyPooledCones = p.Results.identifyPooledCones;
    domainVisualizationLimits = p.Results.domainVisualizationLimits;
    domainVisualizationTicks = p.Results.domainVisualizationTicks;
    visualizedSamplingPositionsGrid = p.Results.visualizedSamplingPositionsGrid;
    visualizedSamplingPositionsGridColor = p.Results.visualizedSamplingPositionsGridColor;
    figureFormat = p.Results.withFigureFormat;

    assert(...
        (size(visualizedSamplingPositionsGridColor,1) == 1) ||...
        (size(visualizedSamplingPositionsGridColor,1) == size(visualizedSamplingPositionsGrid,1)), ...
        'grid colors passed must be either a [1x3], or an [nx3] matrix, where n = # of sampling positions');

    contourGenerationMethod = p.Results.contourGenerationMethod;
	maxNumberOfConesOutsideContour = pStruct.maxNumberOfConesOutsideContour;

	% Load the mosaic
    load(theCenterConnectedMRGCMosaicFullFileName, 'theMRGCMosaic');

    hFig = figure(1000); clf;
    % Prepare figure and axes
    if (isempty(figureFormat))
        ff = PublicationReadyPlotLib.figureComponents('1x1 giant square mosaic');
    else
        ff = figureFormat;
    end
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    theAxes = theAxes{1,1};

    % Compute spatial and chromatic costs
    [theSpatialCompactnessCostsPatch, theSpectralUniformityCostsPatch] = theMRGCMosaic.rfCenterSpatioChromaticCosts();
         
    plotTitle = sprintf('min cone weight visualized: %2.3f',minConeWeightVisualized);
      
    % Visualize the entire mosaic
    if (isempty(domainVisualizationLimits))
        xyMin = min(theMRGCMosaic.rgcRFpositionsDegs,[],1);
        xyMax = max(theMRGCMosaic.rgcRFpositionsDegs,[],1);
        domainVisualizationLimits = [xyMin(1)-0.05 xyMax(1)+0.05 xyMin(2)-0.05 xyMax(2)+0.05];
    else
        fprintf('visualizeFullMosaic: cropping to desired domainVisualizationLimits');
        bandWidth = 0.98*(domainVisualizationLimits(2)-domainVisualizationLimits(1));
        bandHeight = 0.98*(domainVisualizationLimits(4)-domainVisualizationLimits(3));
        bandXo = mean(domainVisualizationLimits(1:2));
        bandYo = mean(domainVisualizationLimits(3:4));
        % Crop it
        theMRGCMosaic.cropToSizeAtEccentricity([bandWidth bandHeight], [bandXo bandYo]);
    end

    if (isempty(domainVisualizationTicks))
        domainVisualizationTicks = struct('x', -50:1:50, 'y', -50:1:50);
    end

    if (~isempty(visualizeRGCsAroundPositionAndWithingPositionDegs))
        visualizedRGCindices = theMRGCMosaic.indicesOfRGCsWithinROI(...
            visualizeRGCsAroundPositionAndWithingPositionDegs(1:2), visualizeRGCsAroundPositionAndWithingPositionDegs(3:4));
    else
        visualizedRGCindices = 1:theMRGCMosaic.rgcsNum;
    end

    

    theMRGCMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', theAxes, ...
        'domainVisualizationLimits', domainVisualizationLimits, ...
        'domainVisualizationTicks', domainVisualizationTicks, ...
        'plottedRFoutlineFaceAlpha', 0.75, ...
        'plottedRFoutlineLineWidth', 1.0, ...
        'plottedRFoutlineFaceColor',  [0 1 0.4], ...
        'minConeWeightVisualized', minConeWeightVisualized, ...
        'maxNumberOfConesOutsideContour', maxNumberOfConesOutsideContour, ...
        'identifyPooledCones', identifyPooledCones, ...
        'identifyInputCones', identifyInputCones, ...
        'pooledConesLineWidth', 1.0, ...
        'visualizedRGCindices', visualizedRGCindices, ...
        'contourGenerationMethod', contourGenerationMethod, ...
        'clearAxesBeforeDrawing', false, ...
        'centerSubregionContourSamples', 20, ...
        'plottedRFoutlineLineWidth', 1.0, ...
        'backgroundColor', [0 0 0], ...
        'plotTitle', plotTitle);

    if (~isempty(visualizedSamplingPositionsGrid))
        hold(theAxes, 'on');
        if (size(visualizedSamplingPositionsGridColor,1) == 1)
            plot(theAxes, visualizedSamplingPositionsGrid(:,1), visualizedSamplingPositionsGrid(:,2), 's', ...
                'MarkerSize', 12, 'MarkerFaceColor', visualizedSamplingPositionsGridColor(1,:)*0.7, ...
                'MarkerEdgeColor', visualizedSamplingPositionsGridColor(1,:), 'LineWidth', 1.0);
        else
            for iPos = 1:size(visualizedSamplingPositionsGrid,1)
                plot(theAxes, visualizedSamplingPositionsGrid(iPos,1), visualizedSamplingPositionsGrid(iPos,2), 's', ...
                    'MarkerSize', 12, 'MarkerFaceColor', visualizedSamplingPositionsGridColor(iPos,:)*0.7, ...
                    'MarkerEdgeColor', visualizedSamplingPositionsGridColor(iPos,:), 'LineWidth', 1.0);
            end
        end
    end

    % Finalize figure using the Publication-Ready format
    ff.box = 'on';
    ff.backgroundColor = [0 0 0];
    PublicationReadyPlotLib.applyFormat(theAxes,ff);

    % Export figure
    theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    thePDFfileName = fullfile(theRawFiguresDir, strrep(theCenterConnectedMRGCMosaicFileName, '.mat', '.pdf' ));
    % Export to PDF
    NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);

    % Also export to PNG for better color saturation
    NicePlot.exportFigToPNG(strrep(thePDFfileName, 'pdf', 'png'),hFig,  300);

end
