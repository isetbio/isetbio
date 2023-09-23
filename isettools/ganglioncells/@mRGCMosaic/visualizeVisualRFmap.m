function visualizeVisualRFmap(theVisualRFmapStruct, retinalRGCRFposDegs, theAxes, varargin)
    % Parse optional input
    p = inputParser;
    p.addParameter('withFigureFormat', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('tickSeparationArcMin', 6, @isscalar);
    p.addParameter('reverseXDir', false, @islogical);
    p.addParameter('gridlessLineWeightingFuncions', false, @islogical);

    p.parse(varargin{:});
    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    reverseXDir = p.Results.reverseXDir;
    gridlessLineWeightingFuncions = p.Results.gridlessLineWeightingFuncions;
    ff = p.Results.withFigureFormat;
    spatialSupportRangeArcMin = tickSeparationArcMin * 4;


     
    % Find the coords of the RF max
    [~,idx] = max(abs(theVisualRFmapStruct.theRFmap(:)));
    [row,col] = ind2sub(size(theVisualRFmapStruct.theRFmap),idx);

    xShiftDegs = mean(theVisualRFmapStruct.spatialSupportDegsX) - theVisualRFmapStruct.spatialSupportDegsX(col);
    yShiftDegs = mean(theVisualRFmapStruct.spatialSupportDegsY) - theVisualRFmapStruct.spatialSupportDegsY(row);
    pixelSizeDegs = theVisualRFmapStruct.spatialSupportDegsX(2)-theVisualRFmapStruct.spatialSupportDegsX(1);
    xShiftSamples = sign(xShiftDegs)*round(abs(xShiftDegs)/pixelSizeDegs);
    yShiftSamples = sign(yShiftDegs)*round(abs(yShiftDegs)/pixelSizeDegs);

    visualRGCRFposDegs(1) = theVisualRFmapStruct.spatialSupportDegsX(col);
    visualRGCRFposDegs(2) = theVisualRFmapStruct.spatialSupportDegsY(row);

    xyLimsArcMin = spatialSupportRangeArcMin/2*[-1 1];
    xyTicksArcMin = 0:(tickSeparationArcMin):60;
    xyTicksArcMin = [-fliplr(xyTicksArcMin(2:end)) xyTicksArcMin];
    xyTicks = 0:(tickSeparationArcMin/60):30;
    xyTicks = [-fliplr(xyTicks(2:end)) xyTicks];

    maxXY = round(spatialSupportRangeArcMin/2);
    spatialSupportDegs = (-maxXY:0.05:maxXY)/60;

    % The cell position
    XLims = visualRGCRFposDegs(1) + [spatialSupportDegs(1) spatialSupportDegs(end)];
    YLims = visualRGCRFposDegs(2) + [spatialSupportDegs(1) spatialSupportDegs(end)];


    surroundBoostedRF = theVisualRFmapStruct.theRFmap;
    surroundBoostedRF(surroundBoostedRF<0) = surroundBoostedRF(surroundBoostedRF<0)*10;

    theRFprofileX = squeeze(sum(theVisualRFmapStruct.theRFmap,1));
    theRFprofileY = squeeze(sum(theVisualRFmapStruct.theRFmap,2));
    maxProfile = max([max(abs(theRFprofileX(:))) max(abs(theRFprofileY(:)))]);

    theRFprofileX = circshift(theRFprofileX ,xShiftSamples);
    theRFprofileY = circshift(theRFprofileY ,yShiftSamples);


    spatialSupportArcMinX = theVisualRFmapStruct.spatialSupportDegsX*60;
    spatialSupportArcMinX = spatialSupportArcMinX - mean(spatialSupportArcMinX);

    spatialSupportArcMinY = theVisualRFmapStruct.spatialSupportDegsY*60;
    spatialSupportArcMinY = spatialSupportArcMinY - mean(spatialSupportArcMinY);


    % Render the RF map
    render2DVisualRF(theAxes{1,1}, ...
        theVisualRFmapStruct.spatialSupportDegsX, theVisualRFmapStruct.spatialSupportDegsY, ...
        theVisualRFmapStruct.theRFmap,  ...
        xyTicks, XLims, YLims, reverseXDir, [], false, ff);

    render1DProfile(theAxes{1,2}, spatialSupportArcMinX, theRFprofileX/maxProfile, ...
        xyTicksArcMin, xyLimsArcMin, reverseXDir, gridlessLineWeightingFuncions, ff);


    render1DProfile(theAxes{1,3}, spatialSupportArcMinY, theRFprofileY/maxProfile, ...
        xyTicksArcMin, xyLimsArcMin, reverseXDir, gridlessLineWeightingFuncions, ff);

    % Render the surround boosted RFmap
    render2DVisualRF(theAxes{1,4}, ...
        theVisualRFmapStruct.spatialSupportDegsX, theVisualRFmapStruct.spatialSupportDegsY, ...
        surroundBoostedRF,  ...
        xyTicks, XLims, YLims, reverseXDir, 'surround boosted (x10)', true, ff);

    drawnow;


    % ============== Export to PLOS directory ==========================
    rawFiguresRoot = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode/Raw';


    % The visual RF map
    pdfFileName = sprintf('InterpolatedRGC_VisualRFmap');

    hFig = figure(100); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    thePLOSAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);
    ax = thePLOSAxes{1,1};

    render2DVisualRF(ax, ...
        theVisualRFmapStruct.spatialSupportDegsX, theVisualRFmapStruct.spatialSupportDegsY, ...
        theVisualRFmapStruct.theRFmap,  ...
        xyTicks, XLims, YLims, reverseXDir, [], false, ff);


    pdfFileNameForPLOS = fullfile(rawFiguresRoot, pdfFileName);
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFig, 300);


    % The visual RF map profile X
    pdfFileName = sprintf('InterpolatedRGC_VisualRFprofileX');

    hFig = figure(200); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    thePLOSAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);
    ax = thePLOSAxes{1,1};

    render1DProfile(ax, spatialSupportArcMinX, theRFprofileX/maxProfile, ...
        xyTicksArcMin, xyLimsArcMin, reverseXDir, gridlessLineWeightingFuncions, ff);

    pdfFileNameForPLOS = fullfile(rawFiguresRoot, pdfFileName);
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFig, 300);


    % The visual RF map profile Y
    pdfFileName = sprintf('InterpolatedRGC_VisualRFprofileY');

    hFig = figure(300); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    thePLOSAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);
    ax = thePLOSAxes{1,1};

    render1DProfile(ax, spatialSupportArcMinY, theRFprofileY/maxProfile, ...
        xyTicksArcMin, xyLimsArcMin, reverseXDir, gridlessLineWeightingFuncions, ff);

    pdfFileNameForPLOS = fullfile(rawFiguresRoot, pdfFileName);
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFig, 300);


end

function render1DProfile(ax, spatialSupportArcMinX, theRFprofile, ...
        xyTicksArcMin, xyLimsArcMin, reverseXDir, gridlessLineWeightingFuncions, ff)

    shadedAreaPlot(ax,spatialSupportArcMinX, theRFprofile, 0, [1 0.85 0.85], [0.7 0 0], 0.5, 1.5);
    hold(ax, 'on');
    plot(ax, spatialSupportArcMinX, theRFprofile*0, 'k-', 'LineWidth', 1.0);

    
    set(ax, 'YLim', [-0.4 1], 'YTick', -1:0.2:1, ...
        'XTick', xyTicksArcMin, 'XLim', xyLimsArcMin, ...
        'XTickLabel', sprintf('%2.1f\n', xyTicksArcMin), ...
        'YTickLabel', {});
    axis(ax, 'square');

    if (gridlessLineWeightingFuncions)
        grid(ax, 'off');
    else
        grid(ax, 'on');
    end
    box(ax, 'on');
    
    xtickangle(ax, 0);

    if (reverseXDir)
        set(ax, 'XDir', 'reverse');
    end

    % Finalize subplot
    if (isempty(ff))
        xlabel(ax, 'space, x (arc min)');
    else
        xlabel(ax, 'space, x (arc min)', 'FontAngle', ff.axisFontAngle);

        % Font size
        set(ax, 'FontSize', ff.fontSize);

        % axis color and width
        set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    end
end

function render2DVisualRF(ax, ...
        spatialSupportDegsX, spatialSupportDegsY, ...
        RFmap,  ...
        xyTicks, XLims, YLims, reverseXDir, plotTitle, noYLabel,  ff)

    rfSensitivityLUT = brewermap(1024, '*RdBu');
    backgroundColor = rfSensitivityLUT(512,:);
    imagesc(ax, spatialSupportDegsX, spatialSupportDegsY, RFmap);

    hold(ax, 'on');
    contour(ax, spatialSupportDegsX, spatialSupportDegsY, RFmap, ...
        exp(-1)*[0.99 1], 'LineWidth', 1.0, 'Color', [0 0 0]);
    contour(ax, spatialSupportDegsX, spatialSupportDegsY, RFmap, ...
        -0.01*[0.99 1], 'LineWidth', 1.0, 'LineStyle', ':', 'Color', [0 0 0]);

    axis(ax, 'image');
    axis(ax, 'xy');
    set(ax, 'Color', backgroundColor, 'CLim', [-1 1], 'XLim', XLims, 'YLim', YLims, ...
        'XTick', xyTicks, 'YTick', xyTicks, ...
        'XTickLabel', sprintf('%2.2f\n', xyTicks), ...
        'YTickLabel', sprintf('%2.2f\n', xyTicks));
    colormap(ax,rfSensitivityLUT);
    xtickangle(ax, 0);
    
    if (reverseXDir)
        set(ax, 'XDir', 'reverse');
    end


    % Finalize subplot
    if (isempty(ff))
        xlabel(ax, 'space, x (deg)');
        if (~noYLabel)
            ylabel(ax, 'space, y (deg)');
        end
    else
        xlabel(ax, 'space, x (deg)', 'FontAngle', ff.axisFontAngle);
        if (~noYLabel)
            ylabel(ax, 'space, y (deg)', 'FontAngle', ff.axisFontAngle);
        end

        if (~isempty(plotTitle))
        title(ax, 'surround boosted (x10)', ...
            'FontSize', ff.titleFontSize, ...
            'FontWeight', ff.titleFontWeight, ...
            'Color', ff.titleColor);
        end

        % Font size
        set(ax, 'FontSize', ff.fontSize);

        % axis color and width
        set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    end
end


function shadedAreaPlot(ax,x,y, baseline, faceColor, edgeColor, faceAlpha, lineWidth)
    x = [x fliplr(x)];
    y = [y y*0+baseline];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, 'FaceAlpha', faceAlpha, 'LineWidth', lineWidth);
end
