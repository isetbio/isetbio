function surroundWeightsMap(pdfExportSubDir, figNo, figPos, ...
	spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
    visualizedSurroundPoolingWeightMapRange, ...
	inputConePositions, inputConeTypes, inputConeWeights, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('axesToRenderIn', [], @(x)(isempty(x)||(isa(x, 'handle'))));
    p.addParameter('domainVisualizationLimits', [], @(x)(isempty(x)||(numel(x)==4)));
    p.addParameter('domainVisualizationTicks', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('noGrid', false, @islogical);

    p.parse(varargin{:});
    axesToRenderIn = p.Results.axesToRenderIn;
    domainVisualizationLimits = p.Results.domainVisualizationLimits;
    domainVisualizationTicks = p.Results.domainVisualizationTicks;

    noGrid = p.Results.noGrid;
    

    spatialSupportRangeArcMin = 3 * spatialSupportTickSeparationArcMin;
    maxXY = round(spatialSupportRangeArcMin/2);
    spatialSupportDegs = (-maxXY:0.05:maxXY)/60;

    dx = (spatialSupportDegs(end)-spatialSupportDegs(1))*0.05;
    XLims = spatialSupportCenterDegs(1) + [spatialSupportDegs(1)-dx spatialSupportDegs(end)+dx];
    YLims = spatialSupportCenterDegs(1) + [spatialSupportDegs(1)-dx spatialSupportDegs(end)+dx];

    if (~isempty(domainVisualizationLimits))
        XLims = domainVisualizationLimits(1:2);
        YLims = domainVisualizationLimits(3:4);
    end

    xTicksDegs = spatialSupportCenterDegs(1) + (-3:3)*spatialSupportTickSeparationArcMin/60;
    yTicksDegs = spatialSupportCenterDegs(2) + (-3:3)*spatialSupportTickSeparationArcMin/60;

    if (~isempty(domainVisualizationTicks))
        xTicksDegs = domainVisualizationTicks.x;
        yTicksDegs = domainVisualizationTicks.y;
    end

    xLabelString = 'ecc, x (degs)';
    yLabelString = 'ecc, y (degs)';

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure no left axis label');

    if (isempty(axesToRenderIn))
        % Initialize figure
        hFig = figure(figNo); clf;
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
        set(hFig, 'Position', [figPos(1) figPos(2) ff.figureSize(1) ff.figureSize(2)]);
        ax = theAxes{1,1};
    else
        ax = axesToRenderIn;
    end

    idx = find(inputConeTypes == cMosaic.LCONE_ID);
    X = inputConePositions(idx,1);
    Y = inputConePositions(idx,2);
    Z = inputConeWeights(idx);
    stem3(ax, X,Y,Z, ...
        'LineWidth', 1.0, ...
        'Color', 'none', ...
        'Marker','o',...
        'MarkerSize', 6, ...
        'MarkerEdgeColor',RGCMosaicConstructor.constants.LcenterColor*0.5,...
        'MarkerFaceColor',RGCMosaicConstructor.constants.LcenterColor);

    hold(ax, 'on');
    idx = find(inputConeTypes == cMosaic.MCONE_ID);
    X = inputConePositions(idx,1);
    Y = inputConePositions(idx,2);
    Z = inputConeWeights(idx);
    stem3(ax, X,Y,Z, ...
        'LineWidth', 1.0, ...
        'Color', 'none', ...
        'Marker','o',...
        'MarkerSize', 6, ...
        'MarkerEdgeColor',RGCMosaicConstructor.constants.McenterColor*0.5,...
        'MarkerFaceColor',RGCMosaicConstructor.constants.McenterColor);


    idx = find(inputConeTypes == cMosaic.SCONE_ID);
    X = inputConePositions(idx,1);
    Y = inputConePositions(idx,2);
    Z = inputConeWeights(idx);
    stem3(ax, X,Y,Z, ...
        'LineWidth', 1.0, ...
        'Color', RGCMosaicConstructor.constants.ScenterColor*0.5, ...
        'Marker','o',...
        'MarkerSize', 6, ...
        'MarkerEdgeColor',RGCMosaicConstructor.constants.ScenterColor*0.5,...
        'MarkerFaceColor',RGCMosaicConstructor.constants.ScenterColor);

    hold(ax, 'off');
    axis(ax, 'square');

    


    if (spatialSupportTickSeparationArcMin/60 > 1-100*eps)
       xTickLabels = sprintf('%2.0f\n', xTicksDegs);
       yTickLabels = sprintf('%2.0f\n', yTicksDegs);
    elseif (spatialSupportTickSeparationArcMin/60>= 0.1-100*eps)
       xTickLabels = sprintf('%2.1f\n', xTicksDegs);
       yTickLabels = sprintf('%2.1f\n', yTicksDegs);
    elseif (spatialSupportTickSeparationArcMin/60 > 0.01-100*eps)
       xTickLabels = sprintf('%2.2f\n', xTicksDegs);
       yTickLabels = sprintf('%2.2f\n', yTicksDegs);
    else
       xTickLabels = sprintf('%2.3f\n', xTicksDegs);
       yTickLabels = sprintf('%2.3f\n', yTicksDegs);
    end

    if (isempty(visualizedSurroundPoolingWeightMapRange))
        visualizedSurroundPoolingWeightMapRange = [0 max(inputConeWeights)];
    end

    set(ax, ...
        'XLim', XLims, 'YLim', YLims, 'ZLim', visualizedSurroundPoolingWeightMapRange, ...
        'XTick', xTicksDegs, 'YTick', yTicksDegs, ...
        'XTickLabel', xTickLabels, ...
        'YTickLabel', yTickLabels, ...
        'ZTickLabel', 0:0.05:1);
   
    legend(ax, {'L-cones', 'M-cones', 'S-cones'}, ...
        'Location', 'NorthOutside', 'Orientation', 'horizontal', 'NumColumns',3);

    xlabel(ax, xLabelString);
    ylabel(ax, yLabelString);

    if (noGrid)
        ff.grid = 'off';
    end

    title(ax, 'surround cone pooling weights map');
    
    elevationAngle = 25;
    azimuthAngle = -45;
    view(ax, azimuthAngle, elevationAngle)

    
    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);

    
    if (isempty(axesToRenderIn))
        % Export figure
        % OLD Way
        %theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
        %thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, pdfFileName);

        pdfFileName = 'surroundWeightsMap.pdf';
        thePDFfileName = fullfile(pdfExportSubDir, pdfFileName);
        NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
    end
end