function visualizeVisualRFmap(theVisualRFmapStruct, retinalRGCRFposDegs, theAxes, varargin)
    % Parse optional input
    p = inputParser;
    p.addParameter('withFigureFormat', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('tickSeparationArcMin', 6, @isscalar);
    p.addParameter('spatialSupportRangeArcMin', [], @isscalar);
    p.parse(varargin{:});
    ff = p.Results.withFigureFormat;
    spatialSupportRangeArcMin = p.Results.spatialSupportRangeArcMin;
    tickSeparationArcMin = p.Results.tickSeparationArcMin;

    if (isempty(spatialSupportRangeArcMin))
        spatialSupportRangeArcMin = 10;
    end

    if (isempty(tickSeparationArcMin))
        tickSeparationArcMin = 2.0;
    end

     
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

    imagesc(theAxes{1,4}, theVisualRFmapStruct.spatialSupportDegsX, theVisualRFmapStruct.spatialSupportDegsY, surroundBoostedRF);
    axis(theAxes{1,4}, 'image');
    axis(theAxes{1,4}, 'xy');
    set(theAxes{1,4}, 'CLim', [-1 1], 'XLim', XLims, 'YLim', YLims, ...
        'XTick', xyTicks, 'YTick', xyTicks, ...
        'XTickLabel', sprintf('%2.2f\n', xyTicks), ...
        'YTickLabel', sprintf('%2.2f\n', xyTicks));
    colormap(theAxes{1,4},brewermap(1024, '*RdBu'));
    xtickangle(theAxes{1,4}, 0);
    

    % Finalize subplot
    if (isempty(ff))
        xlabel(theAxes{1,4}, 'space, x (deg)');
        title(theAxes{1,4}, 'surround boosted x 10');
    else
        xlabel(theAxes{1,4}, 'space, x (deg)', 'FontAngle', ff.axisFontAngle);
        title(theAxes{1,4}, 'surround boosted (x10)', ...
            'FontSize', ff.titleFontSize, ...
            'FontWeight', ff.titleFontWeight, ...
            'Color', ff.titleColor);
        % Font size
        set(theAxes{1,4}, 'FontSize', ff.fontSize);

        % axis color and width
        set(theAxes{1,4}, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    end


    theRFprofileX = squeeze(sum(theVisualRFmapStruct.theRFmap,1));
    theRFprofileY = squeeze(sum(theVisualRFmapStruct.theRFmap,2));
    maxProfile = max([max(abs(theRFprofileX(:))) max(abs(theRFprofileY(:)))]);

    theRFprofileX = circshift(theRFprofileX ,xShiftSamples);
    theRFprofileY = circshift(theRFprofileY ,yShiftSamples);

    spatialSupportArcMinX = theVisualRFmapStruct.spatialSupportDegsX*60;
    spatialSupportArcMinX = spatialSupportArcMinX - mean(spatialSupportArcMinX);

    spatialSupportArcMinY = theVisualRFmapStruct.spatialSupportDegsY*60;
    spatialSupportArcMinY = spatialSupportArcMinY - mean(spatialSupportArcMinY);

    plot(theAxes{1,2}, spatialSupportArcMinX, theRFprofileX/maxProfile, 'r-', 'LineWidth', 1.5);
    set(theAxes{1,2}, 'YLim', [-0.4 1], 'YTick', -1:0.2:1, ...
        'XTick', xyTicksArcMin, 'XLim', xyLimsArcMin, ...
        'XTickLabel', sprintf('%2.1f\n', xyTicksArcMin));
    axis(theAxes{1,2}, 'square');
    grid(theAxes{1,2}, 'on');
    xtickangle(theAxes{1,2}, 0);

    % Finalize subplot
    if (isempty(ff))
        xlabel(theAxes{1,2}, 'space, x (arc min)');
        ylabel(theAxes{1,2}, 'sensitivity');

    else
        xlabel(theAxes{1,2}, 'space, x (arc min)', 'FontAngle', ff.axisFontAngle);
        ylabel(theAxes{1,2}, 'sensitivity', 'FontAngle', ff.axisFontAngle);

        % Font size
        set(theAxes{1,2}, 'FontSize', ff.fontSize);

        % axis color and width
        set(theAxes{1,2}, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    end

    plot(theAxes{1,3}, spatialSupportArcMinY, theRFprofileY/maxProfile, 'r-', 'LineWidth', 1.5);
    set(theAxes{1,3}, 'YLim', [-0.4 1], 'YTick', -1:0.2:1, ...
        'XTick', xyTicksArcMin, 'XLim', xyLimsArcMin, ...
        'XTickLabel', sprintf('%2.1f\n', xyTicksArcMin), ...
        'YTickLabel', {});

    axis(theAxes{1,3}, 'square');
    grid(theAxes{1,3}, 'on');
    xtickangle(theAxes{1,3}, 0);

    % Finalize subplot
    if (isempty(ff))
        xlabel(theAxes{1,3}, 'space, y (arc min)');
    else
        xlabel(theAxes{1,3}, 'space, y (arc min)', 'FontAngle', ff.axisFontAngle);

        % Font size
        set(theAxes{1,3}, 'FontSize', ff.fontSize);

        % axis color and width
        set(theAxes{1,3}, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    end


    
    imagesc(theAxes{1,1}, theVisualRFmapStruct.spatialSupportDegsX, theVisualRFmapStruct.spatialSupportDegsY, theVisualRFmapStruct.theRFmap);
    hold(theAxes{1,1}, 'on');
    contour(theAxes{1,1}, theVisualRFmapStruct.spatialSupportDegsX, theVisualRFmapStruct.spatialSupportDegsY, theVisualRFmapStruct.theRFmap, ...
        exp(-1)*[0.99 1], 'LineWidth', 1.0, 'Color', [0 0 0]);
    contour(theAxes{1,1}, theVisualRFmapStruct.spatialSupportDegsX, theVisualRFmapStruct.spatialSupportDegsY, theVisualRFmapStruct.theRFmap, ...
        -0.01*[0.99 1], 'LineWidth', 1.0, 'LineStyle', '--', 'Color', [0 0 0]);
    axis(theAxes{1,1}, 'image');
    axis(theAxes{1,1}, 'xy');
    colormap(theAxes{1,1},brewermap(1024, '*RdBu'));
    set(theAxes{1,1}, 'CLim', [-1 1], 'XLim', XLims, 'YLim', YLims, ...
        'XTick', xyTicks, 'YTick', xyTicks, ...
        'XTickLabel', sprintf('%2.2f\n', xyTicks), ...
        'YTickLabel', sprintf('%2.2f\n', xyTicks));
    xtickangle(theAxes{1,1}, 0);

    % Finalize subplot
    if (isempty(ff))
        xlabel(theAxes{1,1}, 'space, x (deg)');
        ylabel(theAxes{1,1}, 'space, y (deg)');

    else
        xlabel(theAxes{1,1}, 'space, x (deg)', 'FontAngle', ff.axisFontAngle);
        ylabel(theAxes{1,1}, 'space, y (deg)', 'FontAngle', ff.axisFontAngle);

        % Font size
        set(theAxes{1,1}, 'FontSize', ff.fontSize);

        % axis color and width
        set(theAxes{1,1}, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    end

    drawnow;
end
