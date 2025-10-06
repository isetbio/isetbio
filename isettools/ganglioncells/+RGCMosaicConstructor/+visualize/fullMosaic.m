function [hFig, theAxes] = fullMosaic(theCenterConnectedMRGCMosaicFullFileName, theCenterConnectedMRGCMosaicFileName, ...
        minConeWeightVisualized, pStruct, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('figNo', [], @(x)(isempty(x)||isnumeric(x)));
    p.addParameter('contourGenerationMethod', 'ellipseFitToPooledConePositions', @(x)(ismember(x, mRGCMosaic.validRFsubregionContourGenerationMethods)));
    p.addParameter('visualizeRGCsAroundPositionAndWithingPositionDegs', [], @(x)(isempty(x)||(isnumeric(x)&&(numel(x)==4))));
    p.addParameter('identifyInputCones', false, @islogical);
    p.addParameter('identifyPooledCones', false, @islogical);
    p.addParameter('inputConesAlpha', 0.6, @isscalar);
    p.addParameter('domainVisualizationLimits', [], @(x)(isempty(x)||(isnumeric(x)&&(numel(x)==4))));
    p.addParameter('domainVisualizationTicks', [], @(x)(isempty(x)||isstruct(x)));
    p.addParameter('visualizedSamplingPositionsGrid', [], @(x)(isempty(x)||isnumeric(x)));
    p.addParameter('visualizedSamplingPositionsGridColor', [1 1 0], @(x)(isnumeric(x)&&(size(x,2)==3)));
    p.addParameter('withFigureFormat', [], @isstruct);
    p.addParameter('withoutPlotTitle', false, @islogical);
    p.addParameter('exportVisualizationPDF', true, @islogical);
    p.addParameter('exportVisualizationPNG', true, @islogical);
    p.parse(varargin{:});

    figNo = p.Results.figNo;
    visualizeRGCsAroundPositionAndWithingPositionDegs = p.Results.visualizeRGCsAroundPositionAndWithingPositionDegs;
    identifyInputCones = p.Results.identifyInputCones;
    identifyPooledCones = p.Results.identifyPooledCones;
    inputConesAlpha = p.Results.inputConesAlpha;
    domainVisualizationLimits = p.Results.domainVisualizationLimits;
    domainVisualizationTicks = p.Results.domainVisualizationTicks;
    visualizedSamplingPositionsGrid = p.Results.visualizedSamplingPositionsGrid;
    visualizedSamplingPositionsGridColor = p.Results.visualizedSamplingPositionsGridColor;
    figureFormat = p.Results.withFigureFormat;
    withoutPlotTitle = p.Results.withoutPlotTitle;
    exportVisualizationPDF = p.Results.exportVisualizationPDF;
    exportVisualizationPNG = p.Results.exportVisualizationPNG;

    assert(...
        (size(visualizedSamplingPositionsGridColor,1) == 1) ||...
        (size(visualizedSamplingPositionsGridColor,1) == size(visualizedSamplingPositionsGrid,1)), ...
        'grid colors passed must be either a [1x3], or an [nx3] matrix, where n = # of sampling positions');

    contourGenerationMethod = p.Results.contourGenerationMethod;
	maxNumberOfConesOutsideContour = pStruct.maxNumberOfConesOutsideContour;

	% Load the mosaic
    load(theCenterConnectedMRGCMosaicFullFileName, 'theMRGCMosaic');

    if (isempty(figNo))
        hFig = figure();
    else
        hFig = figure(figNo); 
    end
    clf;
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
         
    if (~withoutPlotTitle)
        plotTitle = sprintf('min cone weight visualized: %2.3f',minConeWeightVisualized);
    else
        plotTitle = '';
    end
      
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
        'inputConesAlpha', inputConesAlpha, ...
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


    if ((exportVisualizationPDF) || (exportVisualizationPNG))
        % Export figure
        theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
        thePDFfileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
            theRawFiguresDir, strrep(theCenterConnectedMRGCMosaicFileName, '.mat', '.pdf' ), ...
            'generateMissingSubDirs', true);
    
        if (exportVisualizationPDF)
            % Export to PDF
            NicePlot.exportFigToPDF(thePDFfileName,hFig,  300, 'beVerbose');
        end

        if (exportVisualizationPNG)
            % Eexport to PNG for better color saturation
            NicePlot.exportFigToPNG(strrep(thePDFfileName, 'pdf', 'png'),hFig,  300, 'beVerbose');
        end
    end

end
